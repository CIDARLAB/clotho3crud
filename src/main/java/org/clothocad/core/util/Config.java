package org.clothocad.core.util;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.Properties;
import lombok.extern.slf4j.Slf4j;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.ParseException;
import org.clothocad.core.ConfigOption;

/** Contains purely static methods */
@Slf4j
public class Config {
    /** Construct CommandLine object from args passed in from main() */
    public static CommandLine parseArgs(String[] args) throws ParseException {
        return new GnuParser().parse(ConfigOption.getOptions(), args);
    }

    /** Construct configuration. Order of precedence is:
      *   1. Command line arguments
      *   2. Configuration file
      *   3. Hard-coded defaults
      *
      * For command line and config file arguments, later arguments take
      * precedence over earlier ones. The default config file name can be
      * overriden on the command line.
      */
    public static Properties get(CommandLine commandLine) {
        final Properties hardCoded = ConfigOption.getDefaultConfig();

        /* Find config file */
        final Properties withoutFile =
            buildCommandLineConfig(commandLine, hardCoded);
        final String configFileName =
            withoutFile.getProperty(ConfigOption.configfile.name());
        if (configFileName == null) {
            log.debug("No config file specified.");
            return withoutFile;
        }

        /* Load config file and build command line arguments */
        return buildCommandLineConfig(
            commandLine,
            loadConfigFile(configFileName, hardCoded)
        );
    }

    /** Return new Properties object from command line arguments and defaults
     */
    private static Properties
    buildCommandLineConfig(CommandLine commandLine,
                           Properties defaults) {
        final Properties out = new Properties(defaults);
        for (final ConfigOption c : ConfigOption.values()) {
            final String key = c.name();
            final String[] values = commandLine.getOptionValues(key);
            if (values != null && values.length > 0)
                out.setProperty(key, values[values.length - 1]);
        }
        return out;
    }

    /** Return new Properties object from config file
     *  If file cannot be read, return defaults.
     */
    private static Properties
    loadConfigFile(String fileName, Properties defaults) {
        final Properties out = new Properties(defaults);
        try (BufferedReader r =
            Files.newBufferedReader(Paths.get(fileName),
                                    StandardCharsets.UTF_8)) {
            out.load(r);
        } catch (IOException e) {
            log.error("Could not load config file at {}: {}", fileName, 
                    e.toString());
            return defaults;
        }
        log.debug("Read configuration file: {}", fileName);
        return out;
    }
}
