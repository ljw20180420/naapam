#/bin/bash

3_to_1() {
    sed 'N;N;s/\n/\t/g'
}

1_to_3() {
    awk -F $'\t' '
        {
            for (i = 1; i < NF - 2; ++i) {
                printf("%s\t", $i)
            }
            for (i = NF - 2; i <= NF; ++i) {
                printf("%s\n", $i)
            }
        }
    '
}

sort_by_count() {
    3_to_1 |
    sort -k2,2gr |
    1_to_3
}

filter_up() {
    3_to_1 |
    awk -F $'\t' -v end1_m_cut1=$1 '
        $8 - $16 == end1_m_cut1{print;}
    ' |
    1_to_3
}

filter_down() {
    3_to_1 |
    awk -F $'\t' -v start2_m_cut2=$1 '
        $11 - $17 == start2_m_cut2{print;}
    ' |
    1_to_3
}

filter_rs() {
    3_to_1 |
    awk -F $'\t' -v rs=$1 '
        $10 == rs{print;}
    ' |
    1_to_3
}
