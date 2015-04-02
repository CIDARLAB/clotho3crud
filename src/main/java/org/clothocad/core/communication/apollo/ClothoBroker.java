package org.clothocad.core.communication.apollo;

import org.apache.activemq.apollo.broker.Broker;
//import org.apache.activemq.apollo.broker.store.leveldb.dto.*;
import org.apache.activemq.apollo.dto.*;
//import org.iq80.leveldb.util.FileUtils;

import java.io.File;
import org.apache.commons.daemon.Daemon;
import org.apache.commons.daemon.DaemonContext;
import org.apache.commons.daemon.DaemonInitException;

public class ClothoBroker 
    implements Daemon {

	private Broker broker = null;
	
	public ClothoBroker() {

        // Creating and initially configuring the broker.
        this.broker = new Broker();
        
        this.broker.setTmp(new File("./tmp"));
        this.broker.setConfig(createConfig());

        this.broker.start(new Runnable() {
            public void run() {
            	System.out.println("The Broker is running...");
            }
        });
        
        Runtime.getRuntime().addShutdownHook(new Thread() {
            public void run() {
            	// delete the data directory
            	//FileUtils.deleteDirectoryContents(new File("./data"));
            }
        });
	}
	
	public void stop() {
        // The broker stops asynchronously. The runnable is invoked once
        // the broker if fully stopped.
        broker.stop(new Runnable(){
            public void run() {
                System.out.println("The broker has now stopped.");
            }
        });
	}
	
    /**
     * Builds a simple configuration model with just plain Java.  Corresponds 1 to 1 with
     * the XML configuration model.  See the Apollo user guide for more details.
     * @return
     */
    private BrokerDTO createConfig() {
        BrokerDTO broker = new BrokerDTO();

        // Brokers support multiple virtual hosts.
        VirtualHostDTO host = new VirtualHostDTO();
        host.id = "localhost";
        host.host_names.add("localhost");
        host.host_names.add("127.0.0.1");
        broker.virtual_hosts.add(host);

        // Control which ports and protocols the broker binds and accepts
        AcceptingConnectorDTO tcpConnector = new AcceptingConnectorDTO();
        tcpConnector.id = "tcp";
        tcpConnector.bind = "tcp://0.0.0.0:61613";
        broker.connectors.add(tcpConnector);

        return broker;
    }

    @Override
    public void init(DaemonContext dc) throws DaemonInitException, Exception {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void start() throws Exception {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void destroy() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
}
