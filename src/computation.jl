include("measures.jl")
include("io.jl")

function lexical_computations(tokens::AbstractArray{String,1})
    l = length(tokens)
    v = length(unique(tokens))

    ttr        = v/l
    guiraud_r  = v/sqrt(l)
    herdan_c   = log(v)/log(l)
    dugast_k   = log(v)/log(log(l))
    dugast_u   = log(l)^2/(log(l) - log(v))
    maas_a2    = (log(l) - log(v))/log(l)^2
    tuldava_ln = (1 - v^2)/(v^2 * log(l))
    brunet_w   = l^(v^-0.172)
    cttr       = v/sqrt(2 * l)
    summer_s   = log(log(v))/log(log(l))

    fs = frequency_spectrum(tokens)
    hapaxlegomena = get(fs, 1, 0)
    dislegomena = get(fs, 2, 0)

    sichel_s = dislegomena/v
    michea_m = v/dislegomena
    honore_h = 100.0 * log(l)/(1 - hapaxlegomena/v)
    ent = sum(freq_size * -log(freq / l) * (freq / l) for (freq, freq_size) in fs)
    yule_k = 10000 * (sum(freq_size * (freq / l)^2 for (freq, freq_size) in fs) - (1 / l))
    simpson_d = sum((freq_size * (freq / l) * ((freq - 1) / (l - 1)) for (freq, freq_size) in fs))
    herdan_vm =  sqrt(sum(freq_size * (freq / l)^2 for (freq, freq_size) in fs) - (1 / v))
    hdd = sum((1 - dhyper(0, freq, l - freq, 42)) / 42 for (word, freq) in fs)
    orlov_z = 1.0 # orlov_z(l, v, fs) # FIXME

    Float64[ttr, guiraud_r, herdan_c, dugast_k, dugast_u, maas_a2,
            tuldava_ln, brunet_w, cttr, summer_s, sichel_s, michea_m, honore_h,
            ent, yule_k, simpson_d, herdan_vm, hdd, orlov_z]
end

"""
Calculates bootstrap for lexical diversity measures as explained in Evert et al. 2017:
http://www.stefan-evert.de/PUB/EvertWankerlNoeth2017.pdf
"""
function bootstrap_computations(tokens::AbstractArray{String,1}, window::Int64)
    N = length(tokens)
    chunks = equal_windows(N, window)
    results = Array{Float64,2}(length(BOOTSTRAP_MEASURE_LIST),length(chunks))
    for (i, start) in enumerate(chunks)
        stop = start + window - 1
        tokens_chunk = @view tokens[start:stop]
        results[:,i] = lexical_computations(tokens_chunk)
    end
    vec(mean(results, 2))
end

function window_computations(tokens::Vector{String}, windows=[100, 1000])
    m = Array{Float64,2}(length(MEASURE_LIST),length(windows))
    mtld_v = mtld(tokens)
    for (i, window) in enumerate(windows)
        m[:,i] = vcat(bootstrap_computations(tokens, window),
                      [mattr(tokens, window); mtld_v; sttr(tokens, window)])
    end
    vec(m)
end

function relabel_measure(measure::String, label::Int64)
    Symbol(string(measure, '_', label))
end

function get_headers(windows::Vector{Int64})
    vcat([Symbol(m) for m in BOOTSTRAP_MEASURE_LIST],
         [relabel_measure(measure, window)
          for window in windows
          for measure in MEASURE_LIST])
end

function compute_all(tokens::Vector{String}, windows=[100, 1000])
    vcat(lexical_computations(tokens), window_computations(tokens, windows))
end

using DataFrames

function compute_corpus!(corpus_dir::String, windows=[100, 1000])
    headers = get_headers(windows)
    files = corpus_files(corpus_dir)
    r = Array{Float64,2}(length(headers),length(files))
    for (i, filename) in enumerate(files)
        x = compute_all(read_tokenized_file(filename), windows)
        r[:,i] = x
    end
    df = DataFrame(vcat(vec(files)'', r'))
    names!(df, hcat("filename", headers))
    df
end
