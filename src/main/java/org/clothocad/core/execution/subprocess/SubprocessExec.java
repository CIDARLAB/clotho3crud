package org.clothocad.core.execution.subprocess;

import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.datums.util.Language;
import org.clothocad.core.util.CloseableProcess;
import org.clothocad.core.util.CloseableRunnable;
import org.clothocad.core.util.CloseableThread;

/** Executes function in an external language interpreter subprocess
 *
 * Our Java process (the "host") communicates with the subprocess over
 * the following channels:
 *   - subprocess standard input:  host --> subprocess messages
 *   - subprocess standard output: subprocess --> host messages
 *   - subprocess standard error:  copied to host; not interpreted
 *
 * Subprocess standard error is captured in a byte array and is not
 * interpreted.
 *
 * Messages are UTF-8 encoded, null-byte terminated JSON values.
 * The types of messages that can be sent to the subprocess are:
 *   function_init:
 *     {"type": "func", "code": <string>, "args": <array>}
 *
 *   api_return:
 *     {"type": "api", "return": <value>}
 *
 *   api_error:
 *     {"type": "api_error", "message": <string>}
 *
 * The messages that can be received from the subprocess are:
 *   function_return:
 *     {"type": "func", "return": <value>}
 *
 *   api_call:
 *     {"type": "api", "name": <string>, "args": <list>}
 *
 * The message protocol and threads of execution are
 * summarized in the following diagram:
 *
 *  Host                    Subprocess
 *main thread
 *  ...
 *   |                       [is born]
 *   |                          |
 *   |                          X (blocks)
 *   |      function_init
 *   | -----------------------> \ (resumes)
 *   X                          |
 *                              |
 *  ...                        ...
 *                              |
 *            api_call          |
 *   / <----------------------- |
 *   |                          X
 *   | api_return or api_error
 *   | -----------------------> \
 *   X                          |
 *                              |
 *  ...                        ...
 *                              |
 *         function_return      |
 *   / <----------------------- |
 *   |                        [dies]
 *  ...
 *
 * The host creates helper threads for reading from the standard output and
 * error streams of the subprocess. When the subprocess dies, so do the helper
 * threads. The diagram above does not show the helper threads.
 */
public class SubprocessExec {
    public static Object
    run(final ServerSideAPI api,
        final Map<String, Object> sourceJSON,
        final List<Object> args,
        final EventHandler eventHandler) {
        
        String slang = (String) sourceJSON.get("language");
        final Language lang = Language.valueOf(slang.toUpperCase());

        try (final CloseableProcess proc =
             ProcessProvider.newProcess(lang);
        ) {
            final Object returnValue;
            final PipeDumper errDumper =
                new PipeDumper(proc.getProcess().getErrorStream());
            final JSONStreamReader reader =
                new JSONStreamReader(proc.getProcess().getInputStream());
            try {
                returnValue = initThreadsAndRun(
                    api,
                    reader,
                    new JSONStreamWriter(proc.getProcess().getOutputStream()),
                    new CloseableRunnable[] {reader, errDumper},
                    (String) sourceJSON.get("code"),
                    args
                );
            } catch (final Exception e) {
                eventHandler.onFail(errDumper.getBytes());
                throw e;
            }
            eventHandler.onSuccess(errDumper.getBytes());
            return returnValue;
        }
    }

    private static Object
    initThreadsAndRun(final ServerSideAPI api,
                      final JSONStreamReader r,
                      final JSONStreamWriter w,
                      final CloseableRunnable[] runnables,
                      final String code,
                      final List<Object> args) {
        final CloseableThread[] threads =
            new CloseableThread[runnables.length];

        for (int i = 0; i < runnables.length; i++)
            threads[i] = new CloseableThread(runnables[i]);

        try (final TryAll tryall = new TryAll(threads)) {
            for (final CloseableThread thread : threads)
                thread.getThread().start();
            return new ExecutionContext(api, r, w).start(code, args);
        }
    }

    public static interface EventHandler {
        void onFail(byte[] standardError);
        void onSuccess(byte[] standardError);
    }
}
