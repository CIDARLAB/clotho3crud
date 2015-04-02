package org.clothocad.core.execution.subprocess;

/* TODO: this whole class is a mess and needs to be re-written.
 * The ServerSideAPI is basically duplicated here.
 */
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.clothocad.core.communication.ServerSideAPI;
import org.clothocad.core.datums.Function;
import org.clothocad.core.datums.ObjectId;

class APIRelayer {
    private static final Map<String, Relayer> RELAYERS = new HashMap<>();
    static {
        RELAYERS.put("get", new GetRelayer());
        RELAYERS.put("set", new SetRelayer());
        RELAYERS.put("run", new RunRelayer());
    }

    static void relay(final ServerSideAPI api,
                      final String name,
                      final List<Object> args,
                      final Callback cb) {
        final Relayer r = RELAYERS.get(name);
        if (r == null)
            throw new RuntimeException(
                String.format("no relayer for '%s'", name)
            );
        r.relay(api, args, cb);
    }

    static interface Callback {
        void onSuccess(Object returnValue);
        void onFail(String message);
    }

    private static interface Relayer {
        void relay(ServerSideAPI api, List<Object> args, Callback cb);
    }

    private static class GetRelayer implements Relayer {
        @Override public void
        relay(final ServerSideAPI api,
              final List<Object> args,
              final Callback cb) {
            if (args.size() != 1)
                throw new RuntimeException();
            cb.onSuccess(api.get(args.get(0)));
        }
    }

    private static class SetRelayer implements Relayer {
        @Override public void
        relay(final ServerSideAPI api,
              final List<Object> args,
              final Callback cb) {
            if (args.size() != 1)
                throw new RuntimeException();
            cb.onSuccess(api.set((Map) args.get(0)));
        }
    }

    private static class RunRelayer implements Relayer {
        @Override public void
        relay(final ServerSideAPI api,
              final List<Object> args,
              final Callback cb) {
            if (args.size() != 2)
                throw new RuntimeException();
            final Object ret;
            final String run_id = (String) args.get(0);
            final List run_args = (List) args.get(1);
            Map obj = new HashMap();
            obj.put("args", run_args);
            obj.put("id", run_id);
            
            try {
                ret = api.run(obj);
            } catch (final Exception e) {
                cb.onFail(e.getMessage());
                return;
            }
            cb.onSuccess(ret);
        }
    }
}
