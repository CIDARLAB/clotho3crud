package org.clothocad.core;

import org.clothocad.core.persistence.DBClassLoader;
import org.clothocad.core.persistence.IdUtils;
import org.clothocad.core.schema.Schema;

import com.google.inject.AbstractModule;
import com.google.inject.name.Names;

import java.util.Properties;

import javax.inject.Singleton;

public class ClothoModule extends AbstractModule {

    protected final Properties config;

    public ClothoModule(Properties config) {
        this.config = config == null ?
            ConfigOption.getDefaultConfig() : config;
    }
    
    public ClothoModule() {
        this(null);
    }

    @Override
    protected void configure() {
        Names.bindProperties(binder(), config);
        requestStaticInjection(Schema.class);
        requestStaticInjection(IdUtils.class);
        bind(DBClassLoader.class).in(Singleton.class);
    }
}
