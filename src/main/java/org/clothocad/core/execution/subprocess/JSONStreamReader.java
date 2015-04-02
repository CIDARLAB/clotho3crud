package org.clothocad.core.execution.subprocess;

import java.io.InputStream;
import java.io.IOException;
import org.clothocad.core.util.ByteArray;
import org.clothocad.core.util.CloseableQueue;
import org.clothocad.core.util.CloseableRunnable;

/** Parses InputStream as JSON values
 *
 * Does blocking reads on an InputStream. Assumes the bytes represent a stream
 * of UTF-8 encoded JSON values. Each JSON value must be terminated with a
 * null byte. Note that the UTF-8 encoding of a JSON value can never contain
 * a null byte.
 *
 * The parsed JSON values are made available via waitValue().
 *
 * Intended for execution on a separate thread.
 */
class JSONStreamReader implements CloseableRunnable {
    private final InputStream pipe;
    private final CloseableQueue parsedValues = new CloseableQueue();
    private final ByteArray buffer = new ByteArray();

    JSONStreamReader(final InputStream pipe) {
        this.pipe = pipe;
    }

    /** Pop one JSON value off queue of read values.
     *
     * Thread-safe.
     * Blocks if queue is empty and stream is still open.
     * Throws NoSuchElementException if queue is empty and stream is dead.
     */
    Object
    waitValue() throws InterruptedException {
        return parsedValues.waitPop();
    }

    @Override public void run() {
        try {
            while (!Thread.interrupted() && loop())
                ;
        } finally {
            parsedValues.close();
        }
    }

    @Override public void
    close() {
        try {
            pipe.close();
        } catch (IOException e) {
        }
    }

    private boolean
    loop() {
        final int b;
        try {
            b = pipe.read();
        } catch (IOException e) {
            return false;
        }
        if (b < 0)
            return false;
        if (b == 0) {
            parsedValues.push(JSONUtil.decodeUTF8(buffer.getArray()));
            buffer.clear();
        } else {
            buffer.add((byte) b);
        }
        return true;
    }
}
