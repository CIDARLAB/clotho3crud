package org.clothocad.core.execution.subprocess;

import java.io.InputStream;
import java.io.IOException;
import org.clothocad.core.util.ByteArray;
import org.clothocad.core.util.CloseableRunnable;

class PipeDumper implements CloseableRunnable {
    private final InputStream pipe;
    private final ByteArray buffer = new ByteArray();

    PipeDumper(final InputStream pipe) {
        this.pipe = pipe;
    }

    byte[] getBytes() {
        return buffer.getArray();
    }

    @Override public void
    run() {
        while (!Thread.interrupted()) {
            final int b;
            try {
                b = pipe.read();
            } catch (IOException e) {
                break;
            }
            if (b < 0)
                break;
            buffer.add((byte) b);
        }
    }

    @Override public void
    close() {
        try {
            pipe.close();
        } catch (IOException e) {
        }
    }
}
