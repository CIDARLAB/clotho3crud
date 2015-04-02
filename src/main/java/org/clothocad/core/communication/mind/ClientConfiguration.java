package org.clothocad.core.communication.mind;

import java.util.Map;
import java.util.HashMap;

/**
 * The visual representation of a client
 */
public class ClientConfiguration {
    void addPage(String socket_id, PageMode mode) {
        Page page = new Page(mode);
        page.setSocketId(socket_id);
        pages.put(socket_id, page);
    }

    void removePage(String socket_id) {
        pages.remove(socket_id);
    }

    public Page getPage(String socket_id) {
        return pages.get(socket_id);
    }

    Iterable<String> getSocketIDs() {
        return pages.keySet();
    }

    /* Map socket_id to Page */
    private final Map<String, Page> pages = new HashMap<String, Page>();
    private Page pageInFocus;
}
