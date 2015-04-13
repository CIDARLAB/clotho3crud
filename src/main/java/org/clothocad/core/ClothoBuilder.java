package org.clothocad.core;

import com.google.inject.Guice;
import com.google.inject.Injector;
import com.google.inject.Key;
import com.google.inject.Module;

import org.apache.shiro.SecurityUtils;
import org.python.google.common.collect.Lists;

import java.util.List;
import java.util.Properties;

public class ClothoBuilder {
    public ClothoBuilder(Properties config, Module... modulesArray) {
        List<Module> modules = Lists.asList(new ClothoModule(config), modulesArray);
        injector = Guice.createInjector(modules);
        SecurityUtils.setSecurityManager(injector.getInstance(org.apache.shiro.mgt.SecurityManager.class));
    }

    public ClothoBuilder(Module... modulesArray) {
        this(null, modulesArray);
    }

    private Injector injector;

    public <T> T get(Class<T> c) {
        return injector.getInstance(c);
    }

    public <T> T get(Key<T> key) {
        return injector.getInstance(key);
    }
}
