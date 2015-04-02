package org.clothocad.core.util;

import java.io.InputStream;

public class CloseableThread implements AutoCloseable {
    private final Thread thread;
    private final AutoCloseable runner;

    public
    CloseableThread(final CloseableRunnable r) {
        runner = r;
        thread = new Thread(r);
    }

    public Thread
    getThread() {return thread;}

    @Override public void
    close() {
        try {
            runner.close();
        } catch (Exception e) {
        }
        thread.interrupt();
        try {
            thread.join();
        } catch (InterruptedException e) {
            /* give up joining */
            thread.interrupt();
        }
    }
}
