def checkHiCParam(paramValue, schema) {
    def jsonSlurper = new groovy.json.JsonSlurper()
    def jsonContent = jsonSlurper.parse ( file ( schema, checkIfExists: true ) )
    def pattern = jsonContent['$defs'].hic_options.properties.hic.pattern
    def match = paramValue ==~ pattern

    return match
}

def checkHiCCombinationsParam(paramValue, schema) {
    def jsonSlurper = new groovy.json.JsonSlurper()
    def jsonContent = jsonSlurper.parse ( file ( schema, checkIfExists: true ) )
    def pattern = jsonContent['$defs'].hic_options.properties.hic_map_combinations.pattern
    def match = paramValue ==~ pattern

    return match
}
