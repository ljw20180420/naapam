rearr() {
    local query=$1
    local ref=$2
    local align_dir="$(dirname "$(dirname "${query}")")/align"
    local align_file="$(basename "${query}")"
    align_file="${align_dir}/${align_file%.query}.alg"
    local line_num="$(wc -l < "${ref}")"
    rearrangement \
        < ${input_file} \
        3< ${ref} |
    gawk -f correct_micro_homology.awk -- \
        ${ref} \
        <(yes up | head -n${line_num}) \
        > ${align_file}
}

query_ref() {
    local query
    local chip
    local ref
    for query in $(find "/home/ljw/sdb1/naapam/query" -name "*.query")
    do
        chip="$(
            basename "${query}" |
            sed -E \
                -e 's/.*/\L&/' \
                -e 's/^.+-(a1|a2|a3|g1n|g2n|g3n)-.+$/\1/'
        )"
        ref="/home/ljw/sdb1/naapam/ref/ref/${chip}.ref"
        printf "%s %s\n" ${query} ${ref}
    done
}

export -f rearr

parallel -a <(query_ref) --jobs 24 rearr
