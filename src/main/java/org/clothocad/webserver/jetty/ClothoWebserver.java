package org.clothocad.webserver.jetty;

import org.clothocad.core.communication.Router;
import org.clothocad.core.communication.ws.ClothoWebSocket;

import com.google.inject.servlet.GuiceFilter;

import lombok.Getter;

import org.eclipse.jetty.security.ConstraintMapping;
import org.eclipse.jetty.security.ConstraintSecurityHandler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.HandlerList;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.server.ssl.SslConnector;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.security.Constraint;
import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.websocket.WebSocketServlet;

import javax.inject.Inject;
import javax.inject.Named;
import javax.servlet.http.HttpServletRequest;

//TODO: convert config to guice module
//TODO: make easy to switch ssl requirement on/off for deploy testing
public class ClothoWebserver {

    @Getter
    final Server server;

    @Inject
    public ClothoWebserver(@Named("port") int nPort,
            @Named("confidentialport") int confidentialPort,
            SslConnector sslConnector,
            @Named("containerServletContext") ServletContextHandler servletHandler,
            final Router router, @Named("clientdirectory") String clientDirectory)
            throws Exception {

        server = new Server();

        //Connectors
        SelectChannelConnector connector0 = new SelectChannelConnector();
        connector0.setPort(nPort);
        connector0.setMaxIdleTime(3600000);
        connector0.setRequestHeaderSize(8192);
        connector0.setConfidentialPort(confidentialPort);
        server.addConnector(connector0);

        sslConnector.setPort(confidentialPort);
        server.addConnector(sslConnector);

        // Connection constraints
        Constraint constraint = new Constraint();
        constraint.setDataConstraint(Constraint.DC_CONFIDENTIAL);

        ConstraintMapping cm = new ConstraintMapping();
        cm.setConstraint(constraint);
        cm.setPathSpec("/*");

        ConstraintSecurityHandler constraintHandler = new ConstraintSecurityHandler();
        constraintHandler.setConstraintMappings(new ConstraintMapping[]{cm});

        // Websocket

        WebSocketServlet wsServlet = new WebSocketServlet() {
            @Override
            public WebSocket doWebSocketConnect(HttpServletRequest request, String protocol) {
                return new ClothoWebSocket(request.getSession().getId(), router); 
            }
        };

        // Static resources
        DefaultServlet staticServlet = new DefaultServlet();

        // Handler stack
        servletHandler.setContextPath("/");
        servletHandler.setResourceBase(clientDirectory);
        servletHandler.setWelcomeFiles(new String[]{"index.html"});

        servletHandler.addFilter(GuiceFilter.class, "/*", null);

        servletHandler.addServlet(new ServletHolder(staticServlet), "/*");
        servletHandler.addServlet(new ServletHolder(wsServlet), "/websocket");
        servletHandler.addServlet(new ServletHolder(new RestApi(router)), "/data/*");

        HandlerList handlers = new HandlerList();
        handlers.addHandler(constraintHandler);
        constraintHandler.setHandler(servletHandler);
        server.setHandler(handlers);
    }

    public void start() throws Exception {
        server.start();
    }
}
