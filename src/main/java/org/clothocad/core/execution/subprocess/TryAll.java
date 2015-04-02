package org.clothocad.core.execution.subprocess;

class TryAll implements AutoCloseable {
    private final AutoCloseable[] closeables;

    TryAll(final AutoCloseable[] closeables) {
        this.closeables = closeables;
    }

    @Override public void
    close() {
        for (final AutoCloseable c : closeables) {
            try {
                c.close();
            } catch (Exception e) {
            }
        }
    }
}
