using ArgParse
using CSV
using DataFrames

include("src/TextComplexityMeasures.jl")
using TextComplexityMeasures

function parse_commandline()
    s = ArgParseSettings()

    @add_arg_table s begin
        "--verbose"
        help = "turn on verbose logging"
        action = :store_true

        "--metadata", "-m"
        help = "TSV file containing metadata for each corpus file"

        "--input-directory", "-i"
        help = "directory containing tokenized files"
        default = "./Tokenized/"
        required = true

        "--out", "-o"
        help = "file where output is saved (leave undefined to print to STDOUT)"
    end

    return parse_args(s)
end

function main()
    parsed_args = parse_commandline()
    if parsed_args["metadata"] !== nothing
        metadata = read_corpus_metadata(parsed_args["metadata"])
    end
    results = compute_corpus!(parsed_args["input-directory"])
    if parsed_args["out"] !== nothing
        CSV.write(parsed_args["out"], results; delim='\t')
    else
        println(results)
    end
end

main()
