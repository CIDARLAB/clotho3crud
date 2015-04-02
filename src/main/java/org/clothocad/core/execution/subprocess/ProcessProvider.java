package org.clothocad.core.execution.subprocess;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.util.CloseableProcess;

class ProcessProvider {
    private static final Map<String, ProcConfig> CONFIGS = new HashMap<>();
    static {
        final Map<String, String> pythonenv = new HashMap<>();
        pythonenv.put("PYTHONPATH",
                      Paths.get("src", "main", "python", "lib").toString());
        CONFIGS.put(
            "Python",
            new ProcConfig(
                new String[] {
                    "python",
                    Paths.get("src", "main", "python", "runner.py").toString()
                },
                pythonenv
            )
        );
    }

    static CloseableProcess newProcess(final Language language) {
        final ProcConfig conf = CONFIGS.get(language.prettyName);
        if (conf.command == null)
            throw new IllegalArgumentException("incorrect language declaration");
        final ProcessBuilder builder = new ProcessBuilder(conf.command);
        builder.environment().putAll(conf.environment);
        try {
            return new CloseableProcess(builder.start());
        } catch (final IOException e) {
            throw new RuntimeException(e);
        }
    }

    private static class ProcConfig {
        private final List<String> command;
        private final Map<String, String> environment;

        private ProcConfig(final String[] cmd, final Map<String, String> env) {
            this.command = Arrays.asList(cmd);
            this.environment = env;
        }
    }
}
