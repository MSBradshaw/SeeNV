// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

def extractSampleInfo(tsvFile) {
    def infos = [] 
    def allLines = tsvFile.readLines()
    for (line in allLines){
        def trimmed = line.trim()
        def cols = trimmed.split()
        
        def idPatient  = cols[0]
        def idSample   = cols[1]
        def sex     = cols[2]
        def bamFile   = returnFile(cols[3])
        def baiFile   = returnFile(cols[4])
        def vcfGzFile   = returnFile(cols[5]) 
        def vcfGzTbiFile   = returnFile(cols[6])
        infos.add([idSample, vcfGzFile, vcfGzTbiFile, bamFile, baiFile])
    }   
    return Channel.from(infos)
}

def singleColumnFileToChannel(tsvFile) {
    def infos = []
    def allLines = tsvFile.readLines()
    for (line in allLines){
        def trimmed = line.trim()
        def cols = trimmed.split()
        infos.add(file(cols[0]))
    }
    return Channel.from(infos)
}

