package org.clothocad.core.communication.mind;

import org.clothocad.core.datums.Doo;

/**
 * When the server tells the client to open a new Page,
 * an AddPageDoo is created to remember this task.
 *
 * @author Kelvin Li
 */
class AddPageDoo extends Doo {
    AddPageDoo(Doo doo, PageMode mode) {
        /* TODO: is this right? */
        super(doo, false);
        this.mode = mode;
    }

    PageMode getPageMode() {
        return mode;
    }

    private PageMode mode;
}
