package org.clothocad.core.execution.subprocess;

import java.io.IOException;
import java.io.OutputStream;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.ServerSideAPI;

/** A helper class for SubprocessExec
 *
 * All this could be easily implemented with purely static methods.
 * We use a class here just for syntactic convenience.
 */
class ExecutionContext {
    private final ServerSideAPI api;
    private final JSONStreamReader reader;
    private final JSONStreamWriter writer;

    ExecutionContext(final ServerSideAPI api,
                     final JSONStreamReader reader,
                     final JSONStreamWriter writer) {
        this.api = api;
        this.reader = reader;
        this.writer = writer;
    }

    /** Communicate with subprocess
     *
     * This method cannot be called more than once.
     * Returns the return value from the function execution
     */
    Object
    start(final String code, final List<Object> args) {
        /* send function_init message */
        final Map<String, Object> value = new HashMap<>();
        value.put("type", "func");
        value.put("code", code);
        value.put("args", args);
        writer.sendValue(value);

        /* do rest of communication */
        return interact();
    }

    /** Interact with subprocess (handle API calls, get return value) */
    private Object
    interact() {
        while (true) {
            final Map value;
            try {
                value = (Map) reader.waitValue();
            } catch (InterruptedException e) {
                throw new RuntimeException(e);
            }
            final String valueType = (String) value.get("type");
            if ("func".equals(valueType)) {
                return handleFunctionReturn(value);
            } else if ("api".equals(valueType)) {
                handleAPICall(value);
            } else {
                /* unrecognized message */
                throw new RuntimeException();
            }
        }
    }

    private Object
    handleFunctionReturn(final Map value) {
        if (!value.containsKey("return"))
            throw new RuntimeException();
        return value.get("return");
    }

    private void
    handleAPICall(final Map value) {
        final Map<String, Object> reply = new HashMap<>();
        final APIRelayer.Callback cb = new APIRelayer.Callback() {
            @Override public void onSuccess(final Object ret) {
                reply.put("type", "api");
                reply.put("return", ret);
            }
            @Override public void onFail(final String message) {
                reply.put("type", "api_error");
                reply.put("message", message);
            }
        };
        APIRelayer.relay(api,
                         (String) value.get("name"),
                         (List) value.get("args"),
                         cb);
        writer.sendValue(reply);
    }
}
