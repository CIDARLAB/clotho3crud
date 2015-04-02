package org.clothocad.core.util;

import java.io.IOException;

public class CloseableProcess implements AutoCloseable {
    private final Process proc;

    public
    CloseableProcess(final Process proc) {
        this.proc = proc;
    }

    public Process
    getProcess() {return proc;}

    @Override public void
    close() {
        try {
            proc.getInputStream().close();
        } catch (IOException e) {
        }
        try {
            proc.getOutputStream().close();
        } catch (IOException e) {
        }
        try {
            proc.getErrorStream().close();
        } catch (IOException e) {
        }
        proc.destroy();
    }
}
