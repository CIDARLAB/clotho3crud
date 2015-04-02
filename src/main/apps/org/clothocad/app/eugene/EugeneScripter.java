package org.clothocad.app.eugene;

import javax.servlet.http.HttpServletRequest;

import org.eclipse.jetty.websocket.WebSocket;
import org.eclipse.jetty.websocket.WebSocketServlet;

public class EugeneScripter extends WebSocketServlet {

	private static final long serialVersionUID = 7273380517601018886L;

	public WebSocket doWebSocketConnect(HttpServletRequest request, String sProtocol) {
		if(sProtocol.equals("eugene")) {
			System.out.println("executing eugene script -> ");
			//return new WSObj();
		} else if(sProtocol.equals("sbol-script")) {
			System.out.println("executing SBOL script -> ");
		} else if(sProtocol.equals("python")) {
			
		} else if(sProtocol.equals("...")) {
			
		}
		return null;
	}
}
