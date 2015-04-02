package org.clothocad.core;

import org.clothocad.core.persistence.Persistor;
import org.clothocad.core.persistence.jongo.JongoModule;
import org.clothocad.core.security.SecurityModule;
import org.clothocad.core.util.JSON;

import com.google.inject.Guice;
import com.google.inject.Injector;

import java.nio.file.Paths;
import java.util.Arrays;
import java.util.Properties;

//Start then navigate to:  http://localhost:8080/#/
public class ClothoStarter extends AbstractClothoStarter {
    public static void main(String[] args) throws Exception {
        baseMain(args, new MainHook() {
            @Override public Injector
            getInjector(Properties config) {
                return Guice.createInjector(
                    new ClothoModule(config),
                    new SecurityModule(),
                    new JongoModule()
                );
            }

            @Override public void call(Injector injector) {
                Persistor persistor = injector.getInstance(Persistor.class);
                ensureMinimalObjects(persistor);
            }
        });
    }

    @Override
    public void start() throws Exception {
        System.out.println("starting with arguments " + Arrays.toString(context.getArguments()));
        main(context.getArguments());
    }

    private static void ensureMinimalObjects(Persistor p) {
        JSON.importTestJSON(
            Paths.get("src", "main", "resources", "json", "essential")
            .toString(), p, false);
    }
}
