using BinDeps
using Compat

@BinDeps.setup

ecos = library_dependency("ecos", aliases=["libecos"])

if is_apple()
    using Homebrew
    provides( Homebrew.HB, "ecos", ecos, os = :Darwin )
end

version = "2.0.5"
win_version = "2.0.2"
provides(Sources, URI("https://github.com/ifa-ethz/ecos/archive/v$version.tar.gz"),
    [ecos], os = :Unix, unpacked_dir="ecos-$version")

prefix = joinpath(BinDeps.depsdir(ecos),"usr")
srcdir = joinpath(BinDeps.depsdir(ecos),"src","ecos-$version")

provides(Binaries, URI("https://cache.julialang.org/https://bintray.com/artifact/download/tkelman/generic/ecos-$win_version.7z"),
    [ecos], unpacked_dir="usr/bin$(Sys.WORD_SIZE)", os = :Windows,
    SHA="b90254220a9a63cba08700f3664519d360f45d363454e5c107e6f30e144a60a1")

# We'll keep this around for emergencies, but OSX users should be able to use Homebrew
provides(SimpleBuild,
    (@build_steps begin
        GetSources(ecos)
        CreateDirectory(joinpath(prefix,"lib"))
        FileRule(joinpath(prefix,"lib","libecos.dylib"),@build_steps begin
            ChangeDirectory(srcdir)
            `make shared`
            `mv libecos.dylib $prefix/lib`
        end)
    end),[ecos], os = :Darwin)

provides(SimpleBuild,
    (@build_steps begin
        GetSources(ecos)
        CreateDirectory(joinpath(prefix,"lib"))
        FileRule(joinpath(prefix,"lib","libecos.so"),@build_steps begin
            ChangeDirectory(srcdir)
            `make shared`
            `mv libecos.so $prefix/lib`
        end)
    end),[ecos], os = :Unix)

@BinDeps.install Dict(:ecos => :ecos)
