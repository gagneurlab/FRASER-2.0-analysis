import pathlib
import wbuild
config['wBuildPath'] =  str(pathlib.Path(wbuild.__file__).parent)

configfile: "wbuild.yaml"
include: config['wBuildPath'] + "/wBuild.snakefile"
include: "Scripts/variant_extraction/variant_extraction.snake"

rule all:
	input: rules.Index.output
