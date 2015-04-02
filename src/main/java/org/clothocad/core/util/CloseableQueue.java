package org.clothocad.core.util;

import java.util.concurrent.BlockingQueue;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.NoSuchElementException;

/** Thread-safe, blocking-pop, closeable queue
 *
 * Handles multiple producers and multiple consumers.
 * The behavior of the queue depends on whether it has been "closed":
 *
 * If queue is "open":
 *   Pop on empty queue blocks until an object is pushed.
 *   Push always succeeds and never waits.
 *
 * If queue is "closed":
 *   Pop on empty queue throws java.util.NoSuchElementException.
 *   Push throws IllegalStateException.
 *
 * The initial state is always "open". Closing is irreversible.
 *
 * This implementation has a push lock but no pop lock.
 */
public class CloseableQueue implements AutoCloseable {
    /* Works around the inability to store null in underyling queue */
    private static final Object NULL_PROXY = new Object();

    /* Tells consumers that queue has been closed */
    private static final Object POISON = new Object();

    /* Tells producers that queue has been closed */
    private boolean is_closed = false;

    /* Underlying queue */
    private final BlockingQueue<Object> queue = new LinkedBlockingQueue<>();

    /** Close this queue
     *
     * This method is idempotent.
     */
    @Override public synchronized void
    close() {
        if (is_closed)
            return;
        is_closed = true;
        queue.add(POISON);
    }

    /** Push object into queue
     *
     * Throws IllegalStateException if queue is closed.
     */
    public synchronized void
    push(final Object obj) {
        if (is_closed)
            throw new IllegalStateException("queue is closed");
        queue.add(obj == null ? NULL_PROXY : obj);
    }

    /** Pop object from queue
     *
     * If queue is empty: blocks if queue is open; throws
     * NoSuchElementException if queue is closed.
     */
    public Object
    waitPop() throws InterruptedException {
        final Object obj = queue.take();
        if (obj != POISON)
            return obj == NULL_PROXY ? null : obj;
        queue.add(POISON);
        throw new NoSuchElementException("queue is closed");
    }
}
