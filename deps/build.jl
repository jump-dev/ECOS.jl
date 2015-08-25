using BinDeps

@BinDeps.setup

ecos = library_dependency("ecos", aliases=["libecos"])

@osx_only begin
    using Homebrew
    provides( Homebrew.HB, "ecos", ecos, os = :Darwin )
end

version = "2.0.2"
provides(Sources, URI("https://github.com/ifa-ethz/ecos/archive/v$version.tar.gz"),
    [ecos], os = :Unix, unpacked_dir="ecos-$version")

prefix = joinpath(BinDeps.depsdir(ecos),"usr")
srcdir = joinpath(BinDeps.depsdir(ecos),"src","ecos-$version")

provides(Binaries, URI("https://cache.e.ip.saba.us/https://bintray.com/artifact/download/tkelman/generic/ecos-$version.7z"),
    [ecos], unpacked_dir="usr$WORD_SIZE/bin", os = :Windows,
    SHA="add47e8b2b14a67c5681a5a77a4dafe0bc4d5efacb38f8c7dffafea79d49d89d")

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

@BinDeps.install [:ecos => :ecos]
