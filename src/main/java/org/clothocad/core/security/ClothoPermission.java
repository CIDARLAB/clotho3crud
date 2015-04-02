/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.security;

import com.google.common.collect.ImmutableSet;
import java.util.Collections;
import java.util.Set;
import static org.clothocad.core.security.ClothoAction.*;

/**
 *
 * @author spaige
 */
public enum ClothoPermission {
    RUN(ImmutableSet.of(run), ImmutableSet.of(run, view, edit, delete, grant)),
    READ(ImmutableSet.of(view, run), ImmutableSet.of(view, edit, delete, grant)),
    WRITE(ImmutableSet.of(view, run, edit), ImmutableSet.of(edit, delete, grant)),
    OWN(ImmutableSet.of(view, run, edit, delete, grant), ImmutableSet.of(delete, grant)),
    PUBLIC(Collections.EMPTY_SET, Collections.EMPTY_SET);
    
    private ClothoPermission(final Set<ClothoAction> actions, Set<ClothoAction> removedActions) {
        this.actions = actions;
        this.removedActions = removedActions;
    }

    public final Set<ClothoAction> actions;
    public final Set<ClothoAction> removedActions;
}
