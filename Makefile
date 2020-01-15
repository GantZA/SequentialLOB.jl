install:
	julia -e 'using Pkg; Pkg.add(PackageSpec(path="$(shell pwd)"))'

compile:
	julia -e 'using PackageCompiler; PackageCompiler.compile_package("SequentialLOB")'

build:
	julia -e 'using PackageCompiler; build_executable("slob_exe.jl", "slobm")'
