/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package org.clothocad.core.util;

import java.io.ByteArrayInputStream;
import java.io.StringWriter;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerException;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;
import org.w3c.dom.Document;
import org.w3c.dom.DocumentType;
import org.w3c.dom.Element;
import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.ErrorHandler;
import org.xml.sax.SAXException;
import org.xml.sax.SAXParseException;

/**
 *
 * @author Admin
 */
public class XMLParser implements ErrorHandler {

    /***************************************************************************/
    /**
     * changes all instances of an attribute in input html to the new value
     * @param html- well formed html formatted as a string
     * @param attribute- an html tag attribute
     * @param prefix-  value that you want to add as the prefix of the current attribute
     * @return html given as a string with the attribute specified set to the new value
     * @throws Exception 
     */
    public static String addPrefixToTagAttribute(String html, String attribute, String prefix) {
        Document xmlDoc = null;
        html = "<html>\n" + html + "\n</html>";
        try {
            DocumentBuilderFactory builderFactory = DocumentBuilderFactory.newInstance();
            builderFactory.setNamespaceAware(true);
            builderFactory.setValidating(true);
            builderFactory.setIgnoringElementContentWhitespace(true);
            DocumentBuilder builder = null;
            builder = builderFactory.newDocumentBuilder();
            builder.setErrorHandler(new XMLParser());
            xmlDoc = builder.parse(new ByteArrayInputStream(html.getBytes()));
            System.out.println("input: " + getStringFromDocument(xmlDoc));
            DocumentType doctype = xmlDoc.getDoctype();
            changeNodeAttribute(xmlDoc.getDocumentElement(), attribute, prefix);
            System.out.println("returning: " + getStringFromDocument(xmlDoc).replaceAll("<html>", "").replaceAll("</html>", ""));
            return getStringFromDocument(xmlDoc).replaceAll("<html>", "").replaceAll("</html>", "");
        } catch (Exception ex) {
//            ex.printStackTrace();
            System.out.println("returning: " + getStringFromDocument(xmlDoc).replaceAll("<html>", "").replaceAll("</html>", ""));
            return getStringFromDocument(xmlDoc).replaceAll("<html>", "").replaceAll("</html>", "");
        }

    }

    static void changeNodeAttribute(Node node, String attribute, String prefix) {
        if (node instanceof Element && node.hasAttributes()) {
            NamedNodeMap attrs = node.getAttributes();
            for (int i = 0; i < attrs.getLength(); i++) {
                ((Element) node).setAttribute(attribute, prefix + ((Element) node).getAttribute(attribute)); //AN ATTRIBUTE IS GETTING CHANGED HERE
            }
        }

        NodeList list = node.getChildNodes();
        if (list.getLength() > 0) {
            for (int i = 0; i < list.getLength(); i++) {
                changeNodeAttribute(list.item(i), attribute, prefix);
            }
        }
    }

    public static String getStringFromDocument(Document doc) {
        try {
            DOMSource domSource = new DOMSource(doc);
            StringWriter writer = new StringWriter();
            StreamResult result = new StreamResult(writer);
            TransformerFactory tf = TransformerFactory.newInstance();
            Transformer transformer = tf.newTransformer();
            transformer.transform(domSource, result);
            return writer.toString();
        } catch (TransformerException ex) {
            ex.printStackTrace();
            return null;
        }
    }

    public void fatalError(SAXParseException spe) throws SAXException {
        System.out.println("Fatal error at line " + spe.getLineNumber());
        System.out.println(spe.getMessage());
        throw spe;
    }

    public void warning(SAXParseException spe) {
        System.out.println("Warning at line " + spe.getLineNumber());
        System.out.println(spe.getMessage());
    }

    public void error(SAXParseException spe) {
        System.out.println("Error at line " + spe.getLineNumber());
        System.out.println(spe.getMessage());
    }
    /***************************************************************************/
}
