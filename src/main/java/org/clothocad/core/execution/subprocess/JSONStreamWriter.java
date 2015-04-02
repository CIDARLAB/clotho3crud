package org.clothocad.core.execution.subprocess;

import java.io.IOException;
import java.io.OutputStream;

class JSONStreamWriter {
    private final OutputStream pipe;

    JSONStreamWriter(final OutputStream pipe) {
        this.pipe = pipe;
    }

    void
    sendValue(final Object value) {
        final byte[] bytes = JSONUtil.encodeUTF8(value);
        try {
            pipe.write(bytes);
            pipe.write(0);
            pipe.flush();
        } catch (IOException e) {
            throw new RuntimeException(e);
        }
    }
}
