package org.clothocad.webserver.jetty;

import org.clothocad.core.communication.Channel;
import org.clothocad.core.communication.Message;
import org.clothocad.core.communication.RestConnection;
import org.clothocad.core.communication.Router;

import org.apache.commons.codec.binary.Base64;
import org.apache.shiro.authz.UnauthorizedException;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.UnsupportedEncodingException;
import java.util.HashMap;
import java.util.Map;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServlet;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

@SuppressWarnings("serial")
public class RestApi extends HttpServlet {

    private static Router router;
    private static Message m, loginMessage, logoutMessage;
    private static Map<String, String> loginMap;
    private static RestConnection rc = new RestConnection("RestConnection");

    // http://stackoverflow.com/questions/15051712/how-to-do-authentication-with-a-rest-api-right-browser-native-clients
    // http://shiro-user.582556.n2.nabble.com/Shiro-and-RESTful-web-services-td5539212.html
    // http://stackoverflow.com/questions/319530/restful-authentication?rq=1
    // http://stackoverflow.com/questions/454355/security-of-rest-authentication-schemes

    // http://shiro.apache.org/
    public RestApi(Router router) {
        this.router = router;
    }

    protected void doGet(HttpServletRequest request, 
        HttpServletResponse response) throws ServletException, IOException {

        // Allows GET calls from domains other than this Clotho's domain
        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");

        if (pathID.length == 0) {
            response.setStatus(HttpServletResponse.SC_OK);
            response.getWriter().write("{\"greeting\": \"Hello Friend!\"}");
            return;
        }

        String[] unamePass = getBasicAuth(request.getHeader("Authorization"));

        login(unamePass);

        String id = pathID[1];

        if (id.equals("query")) {
            String queryString = request.getQueryString();
            Map<String, String> p = splitQuery(queryString);
            m = new Message(Channel.query, p, null, null);
        } else {
            m = new Message(Channel.get, id, null, null);
        }

        try {
            this.router.receiveMessage(this.rc, m);
        } catch (UnauthorizedException ue) {
            response.setStatus(HttpServletResponse.SC_UNAUTHORIZED);
            response.addHeader("WWW-Authenticate", "Basic realm=\"Clotho Rest\"");
            response.addHeader("HTTP/1.0 401", "Unauthorized");
            response.getWriter().write("{\"error\": \"unauthorized access of page\"}");
            return;
        }

        String result = this.rc.getResult().toString();

        logout(unamePass);

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }
    }

    protected void doDelete(HttpServletRequest request, 
        HttpServletResponse response) throws ServletException, IOException {

        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");

        if (pathID.length == 0) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            response.getWriter().write("{\"Required\": \"ID required for delete\"}");
            return;
        }

        String[] unamePass = getBasicAuth(request.getHeader("Authorization"));

        login(unamePass);

        String id = pathID[1];

        m = new Message(Channel.destroy, id, null, null);

        try {
            this.router.receiveMessage(this.rc, m);
        } catch (UnauthorizedException ue) {
            response.setStatus(HttpServletResponse.SC_UNAUTHORIZED);
            response.addHeader("WWW-Authenticate", "Basic realm=\"Clotho Rest\"");
            response.addHeader("HTTP/1.0 401", "Unauthorized");
            response.getWriter().write("{\"error\": \"unauthorized access of page\"}");
            return;
        }

        String result = this.rc.getResult().toString();

        logout(unamePass);

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }
    }

    protected void doPost(HttpServletRequest request, 
        HttpServletResponse response) throws ServletException, IOException {

        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        Map<String, String> p = getRequestBody(request.getReader());

        if (p.isEmpty()) {
            response.getWriter().write("{\"Required\": \"new data to create item with\"}");
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            return;
        }

        String[] unamePass = getBasicAuth(request.getHeader("Authorization"));

        login(unamePass);

        m = new Message(Channel.create, p, null, null);

        try {
            this.router.receiveMessage(this.rc, m);
        } catch (UnauthorizedException ue) {
            response.setStatus(HttpServletResponse.SC_UNAUTHORIZED);
            response.addHeader("WWW-Authenticate", "Basic realm=\"Clotho Rest\"");
            response.addHeader("HTTP/1.0 401", "Unauthorized");
            response.getWriter().write("{\"error\": \"unauthorized access of page\"}");
            return;
        }

        String result = this.rc.getResult().toString();

        logout(unamePass);

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_OK);
        }
    }

    protected void doPut(HttpServletRequest request, 
        HttpServletResponse response) throws ServletException, IOException {

        response.addHeader("Access-Control-Allow-Origin", "*");
        response.setContentType("application/json");

        String[] pathID = request.getPathInfo().split("/");

        if (pathID.length == 0) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            response.getWriter().write("{\"Required\": \"ID required for set\"}");
            return;
        }

        Map<String, String> p = getRequestBody(request.getReader());

        if (p.isEmpty()) {
            response.getWriter().write("{\"Required\": \"new data to set item to\"}");
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
            return;
        }

        String[] unamePass = getBasicAuth(request.getHeader("Authorization"));

        login(unamePass);

        String id = pathID[1];

        if (id.equals("run")) {
            // You can change the url schema, but this makes it a run instead of a set
            m = new Message(Channel.run, p, null, null);
        } else {
            p.put("id", id);
            m = new Message(Channel.set, p, null, null);
        }

        try {
            this.router.receiveMessage(this.rc, m);
        } catch (UnauthorizedException ue) {
            response.setStatus(HttpServletResponse.SC_UNAUTHORIZED);
            response.addHeader("WWW-Authenticate", "Basic realm=\"Clotho Rest\"");
            response.addHeader("HTTP/1.0 401", "Unauthorized");
            response.getWriter().write("{\"error\": \"unauthorized access of page\"}");
            return;
        }

        String result = this.rc.getResult().toString();

        logout(unamePass);

        response.getWriter().write(result);

        if (result.contains("FAILURE")) {
            response.setStatus(HttpServletResponse.SC_BAD_REQUEST);
        } else {
            response.setStatus(HttpServletResponse.SC_CREATED);
        }
    }

    private void login(String[] userPass) {
        if (userPass != null) {
            loginMap = new HashMap<String, String>();
            loginMap.put("username", userPass[0]);
            loginMap.put("password", userPass[1]);
            loginMessage = new Message(Channel.login, loginMap, null, null);
            this.router.receiveMessage(this.rc, loginMessage);
        }
    }

    private void logout(String[] userPass) {
        if (userPass != null) {
            logoutMessage = new Message(Channel.logout, loginMap, null, null);
            this.router.receiveMessage(this.rc, logoutMessage);
        }
    }

    private Map<String, String> getRequestBody(BufferedReader reader) {
        String parts[];
        String keyValue[] = new String[2];
        String currChar = "";
        String key = "";

        Map<String, String> map = new HashMap<String, String>();

        try {
            parts = reader.readLine().split("&");
            for (String kv : parts) {
                keyValue = kv.split("=");
                map.put(keyValue[0], keyValue[1]);
            }
        } catch (IOException ie) {
            map = new HashMap<String, String>();
        } catch (NullPointerException ne) {
            map = new HashMap<String, String>();
        }

        return map;
    }

    private Map<String, String> splitQuery(String queryInfo) {
        Map<String, String> queryPairs = new HashMap<String, String>();
        String[] pairs = queryInfo.split("&");
        for (String pair : pairs) {
            int i = pair.indexOf("=");
            queryPairs.put(pair.substring(0, i), pair.substring(i + 1));
        }
        return queryPairs;
    }

    private String[] getBasicAuth(String authHeader) {
        try {
            String credentials = new String(Base64.decodeBase64(authHeader), "UTF-8");
            return credentials.split(":");
        } catch (UnsupportedEncodingException uee) {
            return null;
        } catch (NullPointerException nee) {
            return null;
        }
    }
}
