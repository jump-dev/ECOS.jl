using BinDeps

@BinDeps.setup

@unix_only begin
    ecos = library_dependency("ecos",aliases=["libecos"])
end

# The last stable version of ECOS seems to be v1.0.3 as of 05/29/2014
# This is safer than unpacking from master which may cause ECOS.jl to
# not work properly
provides(Sources, URI("https://github.com/ifa-ethz/ecos/archive/v1.0.3.zip"),
    [ecos], os = :Unix, unpacked_dir="ecos-1.0.3")

prefix = joinpath(BinDeps.depsdir(ecos),"usr")
srcdir = joinpath(BinDeps.depsdir(ecos),"src","ecos-1.0.3/")

provides(SimpleBuild,
    (@build_steps begin
        GetSources(ecos)
        CreateDirectory(joinpath(prefix,"lib"))
        FileRule(joinpath(prefix,"lib","libecos.dylib"),@build_steps begin
            ChangeDirectory(srcdir)
            `cat ${BinDeps.depsdir(ecos)}/make-dylib.patch` |> `patch Makefile`
            `cat ${BinDeps.depsdir(ecos)}/ecos-fpic.patch` |> `patch ecos.mk`
            `make ecos.dylib`
            `mv libecos.dylib $prefix/lib`
        end)
    end),[ecos], os = :Darwin)

provides(SimpleBuild,
    (@build_steps begin
        GetSources(ecos)
        CreateDirectory(joinpath(prefix,"lib"))
        FileRule(joinpath(prefix,"lib","libecos.so"),@build_steps begin
            ChangeDirectory(srcdir)
            `cat ${BinDeps.depsdir(ecos)}/make-so.patch` |> `patch Makefile`
            `cat ${BinDeps.depsdir(ecos)}/ecos-fpic.patch` |> `patch ecos.mk`
            `make ecos.so`
            `mv libecos.so $prefix/lib`
        end)
    end),[ecos], os = :Unix)

@BinDeps.install [:ecos => :ecos]
