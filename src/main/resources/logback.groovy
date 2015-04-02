//
// Built on Thu Apr 25 21:25:11 CEST 2013 by logback-translator
// For more information on configuration files in Groovy
// please see http://logback.qos.ch/manual/groovy.html

// For assistance related to this tool or configuration files
// in general, please contact the logback user mailing list at
//    http://qos.ch/mailman/listinfo/logback-user

// For professional support please see
//   http://www.qos.ch/shop/products/professionalSupport

def logLevel = TRACE
def rootLevel = INFO

if (System.getProperty("loglevel") != null){
    switch (System.getProperty("loglevel").toUpperCase()){
        case 'OFF':
            logLevel = OFF
            rootLevel = OFF
            break
        case 'ERROR':
            logLevel = ERROR
            rootLevel = ERROR
            break
        case 'WARN':
            logLevel = WARN
            rootLevel = WARN
            break
        case 'INFO':
            logLevel = INFO
            break
        case 'DEBUG':
            logLevel = DEBUG
            break
    }
}

appender("STDOUT", ConsoleAppender) {
  encoder(PatternLayoutEncoder) {
    pattern = "%d{HH:mm:ss.SSS} [%thread] %-5level %logger{36} - %msg%n"
  }
}
root(rootLevel, ["STDOUT"])
logger("org.clothocad", logLevel)
