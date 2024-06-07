library(tidyverse)

taxonomy <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.tax"
fasta <- "benchmarking/trainset19_072023.rdp/trainset19_072023.rdp.fasta"

genera <- read_tsv(taxonomy,
                   col_names = c("accession", "taxonomy")) |>
  mutate(taxonomy = stringi::stri_replace_all_regex(taxonomy, ";$", ""))

fasta_data <- scan(fasta,
                   sep = "\n",
                   what = character(), quiet = TRUE)

sequence_names <- fasta_data[seq(1, length(fasta_data), by = 2)] |>
  str_replace_all(pattern = "^>([^\t]*)\t.*", "\\1")

sequences <- fasta_data[seq(2, length(fasta_data), by = 2)]

seq_table <- tibble(accession = sequence_names, sequence = sequences) |>
  inner_join(genera, by = "accession")

profvis::profvis(
db <- build_kmer_database(seq_table$sequence,
                          seq_table$taxonomy,
                          kmer_size = 8)
)

# from build_kmer_database
sequences = seq_table$sequence
genera = seq_table$taxonomy
kmer_size = 8

genera_indices <- genera_str_to_index(genera)
detected_kmers <- detect_kmers_across_sequences(sequences, kmer_size = kmer_size)
priors <- calc_word_specific_priors(detected_kmers, kmer_size)

# from calc_genus_conditional_prob
detect_list = detected_kmers
genera = genera_indices
word_specific_priors = priors

genus_counts <- tabulate(genera)
n_genera <- length(genus_counts)
n_sequences <- length(genera)
n_kmers <- length(word_specific_priors)
kmer_genus_count <- matrix(0, nrow = n_kmers, ncol = n_genera)
for(i in 1:n_sequences) {
  kmer_genus_count[detect_list[[i]], genera[i]] <-
    kmer_genus_count[detect_list[[i]], genera[i]] + 1
}

bench::mark(
  iterations = 10,
  check = FALSE,
  log((kmer_genus_count + word_specific_priors) / rep(genus_counts + 1, each = n_kmers)),
  calculate_log_probability_cpp(kmer_genus_count, word_specific_priors, genus_counts),
  calculate_log_probability_cpp2(kmer_genus_count, word_specific_priors, genus_counts),
  calculate_log_probability_cpp3(kmer_genus_count, word_specific_priors, genus_counts)
)

# A tibble: 4 × 13
# expression          min median `itr/sec` mem_alloc `gc/sec` n_itr  n_gc total_time result memory
# <bch:expr>       <bch:> <bch:>     <dbl> <bch:byt>    <dbl> <int> <dbl>   <bch:tm> <list> <list>
# 1 log((kmer_gen… 25.95s  29.3s    0.0326    3.79GB  0.0570      4     7      2.05m <NULL> <Rprofmem>
# 2 calculate_log…  22.3s  27.3s    0.0340     1.9GB  0.00849     8     2      3.93m <NULL> <Rprofmem>
# 3 calculate_log…  8.48s  11.9s    0.0816     1.9GB  0.0350      7     3      1.43m <NULL> <Rprofmem>
# 4 calculate_log…  8.64s  11.4s    0.0784     1.9GB  0.0196      8     2       1.7m <NULL> <Rprofmem>
# ℹ 2 more variables: time <list>, gc <list>

unknown_sequence <- sequences[[1]]
bacteroides <- "TACGGAGGATTCAAGCGTTATCCGGATTTATTGGGTTTAAAGGGAGCGTAGGTGGATTGTTAAGTCAGTTGTGAAAGTTTGCGGCTCAACCGTAAAATTGCAGTTGAAACTGGCAGTCTTGAGTACAGTAGAGGTGGGCGGAATTCGTGGTGTAGCGGTTAAATGCTTAGATATCACGAAGAACTCCGATTGCGAAGGCAGCTCACTGGACTGCAACTGACACTGATGCTCGAAAGTGTGGGTATCAAACAGG"
oscillospiraceae <- "TACGTAGGTGGCAAGCGTTGTCCGGATTTACTGGGTGTAAAGGGCGTGTAGCCGGGTTGACAAGTCAGATGTGAAATCCTGCGGCTTAACCGCAGAACTGCATTTGAAACTGTTGATCTTGAGTACTGGAGAGGCAGACGGAATTCCTAGTGTAGCGGTGAAATGCGTAGATATTAGGAGGAACACCAGTGGCGAAGGCGGTCTGCTGGACAGCAACTGACGCTGAAGCACGAAAGTGCGGGGATCGAACAGG"
bacteroidales <- "TACGGAGGATGCGAGCGTTATCCGGATTTATTGGGTTTAAAGGGTGCGTAGGCGGATCGTTAAGTCAGTGGTCAAATTGAGGGGCTCAACCCCTTCCCGCCATTGAAACTGGCGATCTTGAGTGGAAGAGAAGTATGCGGAATGCGTGGTGTAGCGGTGAAATGCATAGATATCACGCAGAACCCCGATTGCGAAGGCAGCATGCCGGCTTCCTACTGACGCTCATGCACGAAAGTGTGGGTAACGAACAGG"

num_bootstraps <- 100
kmer_size <- 8

classify_sequence(unknown = unknown_sequence, database = db,
                  num_bootstraps = num_bootstraps, kmer_size = kmer_size)

filtered <- filter_taxonomy(consensus)
print_taxonomy(filtered)
