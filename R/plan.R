analysis_plan = drake_plan(
  snp_data = read_tsv(
    file = str_glue("{project_folder}/{snp_data_file}"),
    skip = 11) %>%
    clean_names(), # the janitor::clean_names function changes column names so that they are easier to use in R
  
  # Since this is unphased data, there is no difference between a 0/1 and a 1/0 genotype
  # Any subsequent correlation or other analysis is *much* easier if we change the genotypes
  # to integer values, thus
  # 0/0         -> 0
  # 0/1 and 1/0 -> 1
  # 1/1         -> 2
  # I also don't care about the imputation quality values, so I will remove them
  
  snp_data_reformat =
    snp_data %>%
    rename(chrom = number_chrom) %>%
    pivot_longer(cols = contains("_"), # pivoting so that all the RSIDs are in one column and the genotype call is in another
                 names_to = "subject", # simplifies manipulating the values
                 values_to = "call") %>%
    mutate(
      call = str_remove(
        string = call, # To mutate the genotype calls, the quality values first need to be stripped
        pattern = "\\:(.+,?){3}$"
      ) %>% # \\:    -> a literal colon.  colons are special, so we escape it with \\
        # .      -> an alphanumeric or punctuation
        # +      -> match the preceding pattern one or more times
        # ,      -> a literal comma
        # ?      -> match the preceding pattern zero or one times
        # (...)  -> treat everything between the parantheses as a group
        # {3}    -> match the preceding pattern exactly three times
        # $      -> this all needs to be at the end of the string
        # recode takes a character vector and substitues any occurance of the first value with the second
        recode(`0/0` = 0, # so, change 0/0 to 0
               `0/1` = 1, # because of the /, we have to use backticks
               `1/0` = 1, # you can also use a named list, but that list has to have !!! in front of it for weird reasons
               `1/1` = 2
        ),
      subject = str_remove(
        string = subject, # also need to get rid of the things like "SLE|" from the subject values
        pattern = "^[:alpha:]+\\_"
      )
    ),
  
  # We need to make the genotype data wider so that we have a subject-by-snp matrix
  genotypes =
    snp_data_reformat %>%
    select(
      subject,
      id,
      call,
      pos
    ) %>%                           # The only three columns we care about
    filter(pos > lower_bound,
           pos < upper_bound) %>%
    pivot_wider(
      names_from = subject,
      values_from = call
    ) %>%
    column_to_rownames(
      var = "id"
    ) %>%                           # Matrices can only have one data type (in this case, numbers)
    as.matrix(),
  
  # Read the phenotype data (which is an excel file), skipping over the first row
  phenotype = read_xlsx(
    path = str_glue("{project_folder}/{flow_data_file}"),
    skip = 1
  ) %>%
    select(
      subject = `Subject Alias`,                     # Since there are a lot of columns, keep only those of interest
      `Study Year`,
      contains(c("FI", "CD"))                        # The columns with 'FI' and 'CD' appear to be the ones relating to flow data
    ) %>%
    clean_names() %>%
    filter(subject %in% colnames(genotypes)) %>% 
    group_by(subject) %>%                            # We have to deal with the multiple visits.
    top_n(
      n = -1,
      wt = study_year
    ) %>%                                            # I've chosen to do this by only keeping the first visit
    ungroup() %>%
    select(
      -study_year,                                   # We don't need study year any longer
      -blk_median_fi,                                # There were a lot of samples missing values for BLK, CD79A, and the phospho-flow targets
      -percent_of_cd79a,                             # It is either get rid of these analytes or get rid of the subjects
      -matches("^p_")
    ) %>%                                            # because the NAs ruin everything
    mutate_all(as.double) %>%                        # Make sure everything is encoded as an number
    column_to_rownames(var = "subject") %>%
    dplyr::filter(!across(everything(), is.na)) %>%  # NAs are propagated throughout any calculations.  This means that if there are NAs
    as.matrix(),                                     # present when we perform the correlation, the correlation result will be full of NAs.
                                                     # However, that is also true if we calculate the rowSums and we can use that to filter out
                                                     # any rows with NAs
  
  # To perform a correlation between two matrices, we need the number of rows to be equal
  # in this case, we need the subjects to be the same between the genotype and phenotype objects.
  # Also, we can filter out any genotypes where there are only a couple of mutants but making sure
  # the sum of each row is above a certain value.  >5 means that there has to be at least 6 heterozygous
  # or 3 homozygous mutants
  filtered_genotypes = 
    genotypes[
      (rowSums(genotypes) > 5),
      rownames(phenotype)] %>%
    as.data.frame() %>%
    dplyr::filter(!across(everything(), is.na)) %>%
    as.matrix() %>%
    t(),
  
  filtered_phenotype = phenotype %>%
    as.data.frame() %>%
    rownames_to_column(var = 'sample') %>%
    select(sample,
           cd22 = cd22_median_fi,
           cd32 = cd32_median_fi,
           cd72 = cd72_median_fi) %>%
    column_to_rownames(var = 'sample'),
  
  snp_flow_combinations =
    cross2(.x = colnames(filtered_genotypes),
           .y = colnames(filtered_phenotype)) %>%                   # cross2 will make a list of lists that are all possible combinations of snps and flow
    future_map_dfr(.f = as_tibble,                                  # In this case, we end up with 996, 2 member lists (332 * 3)
            .name_repair = ~ set_names(c('snp','flow_analyte'))),   # each of those lists can then be converted into a dataframe/tibble
  
  # we end up with a list of many 2 column/1 row tibbles
  # map_dfr() will take those tibbles and apply bind_rows() to simplify the output
  # (which is what is happening whenever you see sapply(), though in that case, sapply
  # attempts to determine the simplest form automatically, whereas for the map_*()
  # series of functions you have to explicitly state the return type.)
  
  # the first time I attempted to run the pipeline, it essentially hung at the stage of running
  # cor.test.  There are a *lot* of snps in any one chromosome, so I suppose it is not surprising.
  # Instead, divide the possible combinations into 100-line groups and tackle them in parallel
  chunked_df = split(snp_flow_combinations, (seq(nrow(snp_flow_combinations))-1) %/% 100),
  
  
  test_results =
    future_map(chunked_df, function(chunk){
      # map2 lets us walk over two lists simultaneously, so that if we had
      # function f(x, y) and lists A and B, it would apply the
      # function f to each member of A and B like 
      # f(A[[1]], B[[1]]), f(A[[2]], B[[2]]), ... f(A[[n]], B[[n]])
      map2(
        .x = chunk$snp,
        .y = chunk$flow_analyte,
        .f = function(i, j){
          cor.test(
            x = filtered_genotypes[,i],
            y = filtered_phenotype[,j]
          )
        })
    }) %>%
    reduce(c), # merge the individual chunks back together
  
  cor_values = map(test_results, pluck, "estimate") %>% unlist(), # cor.test returns an htest object, we only want two members of those objects
  p_values = map(test_results, pluck, "p.value") %>% unlist(),    # and we want them in a form we can add to a dataframe
  
  # Add the correlation and pvalues to the dataframe with the snp and flow analyte names:
  snp_flow_combinations_with_cor_pval =
    tibble(
      snp = snp_flow_combinations[["snp"]],
      phenotype = snp_flow_combinations[["flow_analyte"]],
      correlation = cor_values,
      pval = p_values
      ) %>%
    replace_na(
      replace = 
        list(
          correlation = 0,
          pval = 0
          )
      ),
  
  # To plot by position, we need to add the position for each snp 
  snp_flow_combinations_with_cor_pval_pos =
    genotypes %>%
    as.data.frame() %>%
    rownames_to_column() %>%
    select(
      rowname,
      pos
    ) %>%
    inner_join(
      snp_flow_combinations_with_cor_pval,
      by = c("rowname" = "snp")
    ) %>%
    rename(snp = rowname),
  
  report = rmarkdown::render(
    input = knitr_in("markdown/report.Rmd"),
    output_file = file_out("reports/report.html"),
    output_dir = "reports",
    quiet = TRUE),
  
  snp_flow_combinations_with_cor_pval_pos %>%
    write_csv(path = file_out("processed_data/phenotype_genotype_correlations.csv"),
              col_names = TRUE)
)
