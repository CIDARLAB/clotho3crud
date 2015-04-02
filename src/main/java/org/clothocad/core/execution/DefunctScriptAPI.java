/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.execution;

import javax.script.ScriptEngine;

/**
 *
 * @author spaige
 */
public interface DefunctScriptAPI {
    public void importFunction(String name);
    public void setEngine(ScriptEngine engine);
}
