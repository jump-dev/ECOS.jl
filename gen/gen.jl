using Clang, ECOS_jll

_GEN_DIR = joinpath(dirname(@__DIR__), "src", "gen")
_COMMON = joinpath(_GEN_DIR, "libecos_common.jl")
function _get_headers(include_path)
    return map(
        f -> joinpath(include_path, f),
        filter(f -> endswith(f, ".h") && f != "timer.h", readdir(include_path)),
    )
end

wc = Clang.init(
    headers = _get_headers(joinpath(ECOS_jll.artifact_dir, "include")),
    output_file = joinpath(_GEN_DIR, "libecos_api.jl"),
    common_file = _COMMON,
    header_wrapped = (root, current) -> root == current,
    header_library = x -> "ecos",
    clang_args = String["-I" * header for header in Clang.find_std_headers()],
    clang_diagnostics = true,
)

run(wc)

function manual_replacements()
    s = read(_COMMON, String)
    for (old, new) in [
        "DBL_MAX + DBL_MAX" => "Inf",
        "ECOS_INFINITY - ECOS_INFINITY" => "NaN",
        "const PRINTTEXT = printf" => "# const PRINTTEXT = printf",
        "const MALLOC = malloc" => "# const MALLOC = malloc",
        "const FREE = free" => "# const FREE = free",
        "const idxint = Cint" => "const idxint = Clong",
    ]
        s = replace(s, old => new)
    end
    write(_COMMON, s)
    return
end

manual_replacements()

rm(joinpath(_GEN_DIR, "LibTemplate.jl"))
