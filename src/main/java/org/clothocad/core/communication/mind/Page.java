package org.clothocad.core.communication.mind;

import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/* A Page is a client browser tab. Every Page has a PageMode,
 * representing multiview, browse mode, trail mode, etc.
 *
 * Page is managed exclusively by ClientConfiguration
 */
public class Page {
    public Page(PageMode mode) {
        this.mode = mode;
    }

    public void addWidget(Widget widget) {
        widgets.put(widget.getId(), widget);
    }

    public void removeWidget(String widget_id) {
        widgets.remove(widget_id);
    }

    boolean toggleVisible(boolean do_toggle) {
        this.is_visible = do_toggle ^ this.is_visible;
        return this.is_visible;
    }

    PageMode getPageMode() {
        return this.mode;
    }

    void setPageMode(PageMode mode) {
        this.mode = mode;
    }

    Widget getWidget(String widget_id) {
        return widgets.get(widget_id);
    }

    Set<String> getWidgetIDs() {
        return widgets.keySet();
    }

    public String getSocketId() {
        return socket_id;
    }

    void setSocketId(String socket_id) {
        this.socket_id = socket_id;
    }
    private boolean is_visible = true;
    private PageMode mode;
    private final Map<String, Widget> widgets = new HashMap<String, Widget>();
    private String socket_id;
    

}
