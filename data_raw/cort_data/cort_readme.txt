# dilution factor for Batch A: 1:100
# for other batches : 1:125


Two subjects each from each batch, from the same pair, separated to be processed in different batch of cort analysis 
(as the plate only fits 37 samples, so each plate has 36 samples = 18 subjects * 2 time points) 

Tried to pick mid-rank animals 


Batch A 
- Pair 4 or Pair 12 ==> Pair 12

Batch B
- Pair 10

Batch D
- Pair 12 

Batch E 
- Pair 10

Batch F 
- 12 

Batch G
- Pair 6 


============================================
$A
# A tibble: 8 x 6
  batch cohort mouseID subjectID pairID ds_rank
  <chr> <chr>  <chr>   <chr>     <chr>    <int>
1 A     A      3       A.3       B            6
2 A     A      4       A.4       W            4
3 A     A      5       A.5       B            7
4 A     A      12      A.12      B            5
5 A     B      4       B.4       B            6
6 A     B      9       B.9       B            7
7 A     B      10      B.10      B            5
8 A     B      12      B.12      W            4

$B
# A tibble: 8 x 6
  batch cohort mouseID subjectID pairID ds_rank
  <chr> <chr>  <chr>   <chr>     <chr>    <int>
1 B     C      6       C.6       B            6
2 B     C      7       C.7       B            4
3 B     C      9       C.9       B            7
4 B     C      10      C.10      W            5
5 B     D      1       D.1       B            7
6 B     D      7       D.7       W            5
7 B     D      10      D.10      B            6
8 B     D      12      D.12      B            4

$D
# A tibble: 8 x 6
  batch cohort mouseID subjectID pairID ds_rank
  <chr> <chr>  <chr>   <chr>     <chr>    <int>
1 D     E      2       E.2       B            4
2 D     E      8       E.8       B            6
3 D     E      9       E.9       W            7
4 D     E      12      E.12      W            5
5 D     F      7       F.7       B            6
6 D     F      8       F.8       W            5
7 D     F      10      F.10      B            4
8 D     F      12      F.12      B            7

$E
# A tibble: 8 x 6
  batch cohort mouseID subjectID pairID ds_rank
  <chr> <chr>  <chr>   <chr>     <chr>    <int>
1 E     G      6       G.6       W            4
2 E     G      8       G.8       B            7
3 E     G      9       G.9       W            6
4 E     G      10      G.10      B            5
5 E     H      2       H.2       W            5
6 E     H      5       H.5       B            4
7 E     H      7       H.7       W            6
8 E     H      10      H.10      W            7

$F
# A tibble: 8 x 6
  batch cohort mouseID subjectID pairID ds_rank
  <chr> <chr>  <chr>   <chr>     <chr>    <int>
1 F     I      3       I.3       B            6
2 F     I      6       I.6       B            5
3 F     I      8       I.8       B            4
4 F     I      12      I.12      B            7
5 F     J      1       J.1       W            6
6 F     J      3       J.3       W            5
7 F     J      4       J.4       B            7
8 F     J      12      J.12      W            4

$G
# A tibble: 8 x 6
  batch cohort mouseID subjectID pairID ds_rank
  <chr> <chr>  <chr>   <chr>     <chr>    <int>
1 G     K      2       K.2       W            6
2 G     K      4       K.4       B            4
3 G     K      6       K.6       W            5
4 G     K      1       K.1       W            7
5 G     L      4       L.4       W            4
6 G     L      6       L.6       B            6
7 G     L      8       L.8       B            7
8 G     L      5       L.5       W            5
