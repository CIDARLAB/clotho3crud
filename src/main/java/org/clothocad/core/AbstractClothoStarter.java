package org.clothocad.core;

import org.clothocad.core.util.Config;
import org.clothocad.webserver.jetty.ClothoWebserver;

import com.google.inject.Injector;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.ParseException;
import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;

import java.util.Properties;

abstract public class AbstractClothoStarter implements Daemon {
    /** This is part of a mechanism to factor out common code
      * across multiple main classes.
      */
    protected static interface MainHook {
        Injector getInjector(Properties config);
        void call(Injector injector);
    }

    private static ClothoWebserver server;
    protected DaemonContext context;

    /** Generic program entry point.
      * Customized setup is passed in the hook parameter.
      */
    protected static void baseMain(String[] args, MainHook hook) throws Exception {
        final CommandLine commandLine;
        try {
            commandLine = Config.parseArgs(args);
        } catch (ParseException e) {
            System.out.println(e.getMessage());
            printHelp();
            return;
        }
        if (commandLine.hasOption("help")) {
            printHelp();
            return;
        }
        final Properties config = Config.get(commandLine);
        /* Do custom setup */
        /* TODO: if keystorepass option passed without arg,
         * then prompt for password
         */
        Injector injector = hook.getInjector(config);
        hook.call(injector);

        server = injector.getInstance(ClothoWebserver.class);
        server.start();
    }

    private static void printHelp() {
        new HelpFormatter().printHelp(" ", ConfigOption.getOptions());
    }

    @Override
    public void init(DaemonContext dc) {
        context = dc;
    }

    @Override
    public void stop() throws Exception {
        server.getServer().stop();
    }

    @Override
    public void destroy() {
        System.out.println("done.");
    }
}
