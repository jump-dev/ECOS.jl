using BinDeps

@BinDeps.setup

@unix_only begin
    ecos = library_dependency("ecos",aliases=["libecos"])
end

provides(Sources, URI("https://github.com/ifa-ethz/ecos/archive/master.zip"),
    [ecos], os = :Unix, unpacked_dir="ecos-master")

prefix = joinpath(BinDeps.depsdir(ecos),"usr")
srcdir = joinpath(BinDeps.depsdir(ecos),"src","ecos-master/")

provides(SimpleBuild,
    (@build_steps begin
        GetSources(ecos)
        CreateDirectory(joinpath(prefix,"lib"))
        FileRule(joinpath(prefix,"lib","libecos.so"),@build_steps begin
            ChangeDirectory(srcdir)
            `cat ${BinDeps.depsdir(ecos)}/make-so.patch` |> `patch Makefile`
            `make CFLAGS="-fPIC" ecos.so`
            `mv libecos.so $prefix/lib`
        end)
    end),[ecos], os = :Unix)

@BinDeps.install [:ecos => :ecos]
