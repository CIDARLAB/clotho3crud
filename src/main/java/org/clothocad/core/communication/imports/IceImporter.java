/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.communication.imports;

import com.fasterxml.jackson.core.JsonParser;
import com.fasterxml.jackson.databind.ObjectMapper;
import com.fasterxml.jackson.xml.XmlMapper;
import com.google.common.base.Joiner;
import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.OutputStreamWriter;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import javax.net.ssl.HttpsURLConnection;
import org.clothocad.core.persistence.Persistor;
import org.clothocad.model.BasicPart;
import org.clothocad.model.FreeForm;
import org.clothocad.model.Part;
import org.clothocad.model.Person;

/**
 *
 * @author spaige
 */
public class IceImporter {

    public static String httpSend(HttpURLConnection conn, String requestBody) throws IOException {
        conn.setDoOutput(true);
        conn.setDoInput(true);
        OutputStreamWriter writer = new OutputStreamWriter(conn.getOutputStream());
        writer.write(requestBody);
        writer.close();
        if (conn.getResponseCode() == HttpURLConnection.HTTP_OK) {
            return new BufferedReader(new InputStreamReader(conn.getInputStream())).readLine();
        } else {
            throw new RuntimeException("got response code " + new Integer(conn.getResponseCode()).toString());
        }

    }

    public static String httpGet(HttpURLConnection conn) throws IOException {
        conn.setDoInput(true);
        StringBuilder out = new StringBuilder();
        String aux;
        BufferedReader reader = new BufferedReader(new InputStreamReader(conn.getInputStream()));
        while ((aux = reader.readLine()) != null) {
            out.append(aux);
        }
        return out.toString();
    }

    public static List<BasicPart> importData(List<Integer> entries) {
        List<BasicPart> out = new ArrayList<>();
        try {
            //login
            URL ice = new URL("https://public-registry.jbei.org/gwt_ice/ice");
            HttpsURLConnection client = (HttpsURLConnection) ice.openConnection();
            client.setRequestMethod("POST");
            client.setRequestProperty("Content-Type", "text/x-gwt-rpc; charset=UTF-8");
            client.setRequestProperty("X-GWT-Module-Base", "https://public-registry.jbei.org/gwt_ice/");
            client.setRequestProperty("X-GWT-Permutation", "13BBE744E36634F5EBE2ABA886E06E57");
            String response = httpSend(client, "7|0|7|https://public-registry.jbei.org/gwt_ice/|AA9EDD4C8A99C9E2FBB951C601A209D2|org.jbei.ice.client.RegistryService|login|java.lang.String/2004016611|stmpaige@gmail.com|blue j 5|1|2|3|4|2|5|5|6|7|");
            ObjectMapper mapper = new ObjectMapper();
            mapper.configure(JsonParser.Feature.ALLOW_SINGLE_QUOTES, true);
            List r = mapper.readValue(response.substring(4), List.class);

            String sessionid = (String) ((List) r.get(16)).get(8);

            //get exported entries
            URL queryURL = new URL("https://public-registry.jbei.org/export?type=XML&entries=" + Joiner.on(",").join(entries));
            client = (HttpsURLConnection) queryURL.openConnection();
            client.setRequestMethod("GET");
            client.setRequestProperty("Cookie", "gd-ice=" + sessionid);
            response = httpGet(client);

            XmlMapper xmlMapper = new XmlMapper();
            List<HashMap> entities = xmlMapper.readValue(response, List.class);

            //parse and save
            for (HashMap entity : entities) {
                if (entity.containsKey("recordType")) {
                    out.add(parseJBEIEntity(entity));
                }
            }

        } catch (MalformedURLException ex) {
            ex.printStackTrace();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
        return out;
    }

    public static void directImport(Persistor p, List<Integer> entries, int attempt) {
        for (BasicPart entry : importData(entries)) {
            p.save(entry);
        }
    }

    private static BasicPart parseToPart(HashMap entity) {
        //String name, String shortdescription, String seq, Format form, Person author
        String name;
        String seq;
        Person author;
        String shortdescription;

        seq = parseSeq(entity);
       

        author = parseAuthor(entity);

        //XXX: 
        name = parsePartName(entity);

        shortdescription = parseDescription(entity);
        
        BasicPart part = new BasicPart(name, shortdescription, seq, new FreeForm(), author);
        
        return part;
    }

    private static String parseSeq(HashMap entity) {
        try {
            return ((HashMap) entity.get("seq")).get("sequence").toString();
        } catch (NullPointerException e) {
            return null;
        }
    }

    private static String parsePartName(HashMap entity) {
        try {
            return ((HashMap) entity.get("partNames")).get("partName").toString();
        } catch (NullPointerException e) {
            return null;
        }
    }

    private static Person parseAuthor(HashMap entity) {
        try {
            return new Person(((HashMap) entity.get("creator")).get("personName").toString());
        } catch (NullPointerException e) {
            return null;
        }
    }

    private static String parseDescription(HashMap entity) {
        try {
            return entity.get("longDescription").toString();
        } catch (NullPointerException e) {
            return null;
        }
    }

    private static BasicPart parseJBEIEntity(HashMap entity) {
        switch (entity.get("recordType").toString()) {
            case "part":
                return parseToPart(entity);
            default:
                return null;
        }
    }
}
