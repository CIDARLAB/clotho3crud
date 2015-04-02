package org.clothocad.core.communication.mind;


public class Widget {
//    public Widget (Page page, View view) {
//        this.widget_id = java.util.UUID.randomUUID().toString();
//        this.position.page = page;
//        this.view = view;
//    }

    public String getId () {
        return this.widget_id;
    }

    void setPositionAbsolute (Double x, Double y) {
        this.position.sibling = null;
        this.position.x = x;
        this.position.y = y;
    }

    void setDimensions (Double w, SizeType wt, Double h, SizeType ht) {
        this.size.width = w;
        this.size.width_type = wt;
        this.size.height = h;
        this.size.height_type = ht;
    }
//
//    public View getView() {
//        return view;
//    }

    public Page getPage() {
        return position.page;
    }

    static enum SiblingRelation {
        ABOVE,
        BELOW,
    }

    static enum SizeType {
        RELATIVE,
        ABSOLUTE,
    }

    private static String getWidgetRootID () {
        return "widget_space_root";
    }

    private class Position {
        Page page;
        Widget parent;

        /* for absolute position */
        Double x = null;
        Double y = null;

        /* for sibling-relative position */
        Widget sibling = null;

        /* relation of self _to_ sibling */
        SiblingRelation relation = SiblingRelation.BELOW;
    }

    private class Size {
        Double width = null;
        SizeType width_type = SizeType.RELATIVE;
        Double height = null;
        SizeType height_type = SizeType.RELATIVE;
    }
    

//    private final View view;
    private String widget_id;
    private Position position = new Position();  //Include Page ref
    private Size size = new Size();

}
