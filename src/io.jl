"""
File input/output.
"""

read_corpus_metadata(filename) = CSV.read(filename; delim='\t')

read_tokenized_file(filename) = String[
    token
    for line in eachline(filename)
    for token in split(line, " ")
    if token != "" && token != "<EOS>" && token != "<PGB>"
]

corpus_files(dir) = String[
    joinpath(root, file)
    for (root, _, files) in walkdir(expanduser(dir); follow_symlinks=true)
    for file in files
    if splitext(file)[2] == ".txt"
]

read_corpus(files) = map(read_tokenized_file, files)
