using DataFrames
using DataArrays
import Rmath: dhyper

"""
Lexical complexity measures.

These and other measures can be found at http://pers-www.wlv.ac.uk/~in4326/papers/\$U50.pdf
Note: brunet_w is not specified correctly there (missing nested ^)

Sample size and vocabulary size:
l = text_length
v = vocabulary_size
"""

ttr(l, v)        = v/l
guiraud_r(l, v)  = v/sqrt(l)
herdan_c(l, v)   = log(v)/log(l)
dugast_k(l, v)   = log(v)/log(log(l))
dugast_u(l, v)   = log(l)^2/(log(l) - log(v))
maas_a2(l, v)    = (log(l) - log(v))/log(l)^2
tuldava_ln(l, v) = (1 - v^2)/(v^2 * log(l))
brunet_w(l, v)   = l^(v^-0.172)
cttr(l, v)       = v/sqrt(2 * l)
summer_s(l, v)   = log(log(v))/log(log(l))

"""
Fequency spectrum and vocabulary size

fs = frequency_spectrum
l = text_length
v = vocabulary_size

Notes:
-   fs[1] = hapax legomena (number of tokens that occur once)
-   fs[2] = dislegomena (number of tokens that occur twice)
"""

function frequency_spectrum(tokens)
    fl = Dict{String,Int64}()
    fs = Dict{Int64,Int64}()
    for token in tokens
        fl[token] = get(fl, token, 0) + 1
    end
    for f in values(fl)
        fs[f] = get(fs, f, 0) + 1
    end
    fs
end

sichel_s(v, fs) = get(fs, 2, 0)/v
michea_m(v, fs) = v/get(fs, 2, 0)
honore_h(l, v, fs) = 100.0 * log(l)/(1 - get(fs, 1, 0)/v)
ent(l, fs) = sum(freq_size * -log(freq / l) * (freq / l) for (freq, freq_size) in fs)
yule_k(l, fs) = 10000 * (sum(freq_size * (freq / l)^2 for (freq, freq_size) in fs) - (1 / l))
simpson_d(l, fs) = sum((freq_size * (freq / l) * ((freq - 1) / (l - 1)) for (freq, freq_size) in fs))
herdan_vm(l, v, fs) =  sqrt(sum(freq_size * (freq / l)^2 for (freq, freq_size) in fs) - (1 / v))
hdd(l, fs, sample_size=42) = sum((1 - dhyper(0, freq, l - freq, sample_size)) / sample_size for (word, freq) in fs)

"""
Parameters of probabilistic models
"""

function orlov_z(l, v, fs, max_iterations=100, min_tolerance=1.0)
    most_frequent = max(keys(fs)...)
    p_star = most_frequent / l
    lower_z, upper_z = nothing, nothing
    z = fld(l, 100)  # our initial guess
    found = false
    for i in range(1, max_iterations)
        estimated_v = (z / log(p_star * z)) * (l / (l - z)) *log(l / z)
        if abs(v - estimated_v) <= min_tolerance
            found = true
            break
        end
        if estimated_v < v
            lower_z = z
            if upper_z != nothing
                z = fld(z + upper_z, 2)
            else
                z *= 2
            end
        else
            upper_z = z
            if lower_z != nothing
                z = fld(z + lower_z, 2)
            else
                z = fld(z,  2)
            end
        end
    end
    if !found
        println("Exceeded max_iterations")
    end
    z
end


# Measure using whole text


"""
Calculate the Moving-Average Type-Token Ratio (Covington and McFall, 2010).

Source:
M.A. Covington, J.D. McFall: Cutting the Gordon Knot. In: Journal of Quantitative Linguistics 17,2 (2010), p. 94-100. DOI: 10.1080/09296171003643098
"""

function mattr(tokens, window=1000)
    ttr_values = Vector{Float64}()
    window_frequencies = Dict{String, Int64}()
    for token in tokens
        window_frequencies[token] = get(window_frequencies, token, 0) + 1
    end
    for window_start in 1:length(tokens) - (window + 2)
        window_end = window_start + window
        word_to_pop = tokens[window_start]
        window_frequencies[word_to_pop] -= 1
        window_frequencies[tokens[window_end]] += 1
        if window_frequencies[word_to_pop] == 0
            delete!(window_frequencies, word_to_pop)
        end
        push!(ttr_values, length(window_frequencies) / window)
    end
    mean(ttr_values)
end


"""
MTLD is an index of a text’s LD, evaluated sequentially. It is calculated as the mean length of sequential word strings in a text that maintain a given TTR value (here, .720).

Implementation following the description in McCarthy and Jarvis (2010): https://link.springer.com/article/10.3758/BRM.42.2.381
"""
function mtld_calculation(tokens, factor_size)
    factors = 0
    #factor_lengths = Vector{Int}()
    types = Set{String}()
    token_count = 0
    for token in tokens
        push!(types, token)
        token_count += 1
        ttr = length(types) / token_count
        if ttr <= factor_size
            factors += 1
            #push!(factor_lengths, token_count)
            types = Set{String}()
            token_count = 0
        end
    end
    if token_count > 0
        ttr = length(types) / token_count
        factors += (1 - ttr) / (1 - factor_size)
        #push!(factor_lengths, token_count)
    end
    return length(tokens) / factors
end

function mtld(tokens, factor_size=0.72)
    mean((mtld_calculation(tokens, factor_size),
          mtld_calculation(reverse(tokens), factor_size)))
end

"""
Standardized Type-Token Ratio

The last partial chunk, if present, is ignored.
"""

function equal_windows(N, window)
    leftover = mod(N, window)
    if leftover == 0
        1:window:N
    else
        1:window:N-leftover
    end
end

function sttr(tokens, window)
    N = length(tokens)
    chunks = equal_windows(N, window)
    sttr_chunks = Vector{Float64}(length(chunks))
    for start in chunks
        stop = start + window - 1
        tokens_chunk = tokens[start:stop]
        ttr = length(unique(tokens_chunk)) / length(tokens_chunk)
        push!(sttr_chunks, ttr)
    end
    ci = 1.96 * std(sttr_chunks) / sqrt(length(sttr_chunks))
    mean(sttr_chunks) #, ci # TODO how to deal with ci's for bootstrapped measures too
end

const BOOTSTRAP_MEASURE_LIST = String[
    "ttr",
    "guiraud_r",
    "herdan_c",
    "dugast_k",
    "dugast_u",
    "maas_a2",
    "tuldava_ln",
    "brunet_w",
    "cttr",

    "summer_s",
    "sichel_s",
    "michea_m",
    "honore_h",
    "ent",
    "yule_k",
    "simpson_d",
    "herdan_vm",
    "hdd",

    "orlov_z",
]

const NON_BOOTSTRAP_MEASURE_LIST = String[
    "mattr",
    "mtld",
    "sttr",
]

const MEASURE_LIST = vcat(BOOTSTRAP_MEASURE_LIST, NON_BOOTSTRAP_MEASURE_LIST)
