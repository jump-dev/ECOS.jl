using BinDeps

@BinDeps.setup

ecos = library_dependency("ecos", aliases=["libecos"])

#@osx_only begin
#    using Homebrew
#    provides( Homebrew.HB, "ecos", ecos, os = :Darwin )
#end

# This is the git commit that includes our merged patches as of 07/30/2014
# This is safer than unpacking from master which may cause ECOS.jl to
# not work properly
provides(Sources, URI("https://github.com/ifa-ethz/ecos/archive/d206a556a83396756f3200964de162b4a7523c62.tar.gz"),
    [ecos], os = :Unix, unpacked_dir="ecos-d206a556a83396756f3200964de162b4a7523c62")

prefix = joinpath(BinDeps.depsdir(ecos),"usr")
srcdir = joinpath(BinDeps.depsdir(ecos),"src","ecos-d206a556a83396756f3200964de162b4a7523c62")

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
