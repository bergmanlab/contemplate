
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.w3c.dom.NamedNodeMap;

import org.xml.sax.SAXException;
import org.xml.sax.ContentHandler;
//import org.xml.sax.Attributes;
import org.xml.sax.helpers.AttributesImpl;
import org.xml.sax.SAXNotRecognizedException;
/**
 * A minimal class that traverses a part of a DOM tree
 * and generates SAX events from it to allow it to be handled
 * by a SAX ContentHandler.
 *
 * @author David Huen
 */
public class DOM2SAXConverter
{
    private static final String NAMESPACESURI = "http://xml.org/sax/features/namespaces";
    private static final String NAMESPACEPREFIXURI = "http://xml.org/sax/features/namespace-prefixes";

    boolean supplyLocalName = true;
    boolean supplyPrefix = false;

    /**
     * Process the specified subtree and generate
     * events from it.
     */
    static public void convertToSAX(Node root, ContentHandler handler)
        throws SAXException
    {
        switch (root.getNodeType()) {
            case Node.ELEMENT_NODE:
                handleElement(root, handler);
                break;

            case Node.TEXT_NODE:
            case Node.CDATA_SECTION_NODE:
                handleText(root, handler);
                break;

            case Node.DOCUMENT_NODE:
                break;

            case Node.PROCESSING_INSTRUCTION_NODE:
                break;

            case Node.ENTITY_REFERENCE_NODE:
                break;

            case Node.DOCUMENT_TYPE_NODE:
                break;
        }
    }

    static private void handleElement(Node node, ContentHandler handler)
        throws SAXException
    {
        //  get the namespace URI
        String nsURI = node.getNamespaceURI();
        String localName = node.getNodeName();

        // construct the FQN
        String qName;
/*        if (nsURI != null) 
            qName = nsURI + localName;
        else 
            nsURI = qName = null;
*/
        if (nsURI == null) nsURI = "";
        qName = nsURI + localName;

        // construct attributes
        AttributesImpl attr = new AttributesImpl();

        NamedNodeMap attrMap = node.getAttributes();

        for (int i=0; i < attrMap.getLength(); i++) {
            Node att = attrMap.item(i);

            String attNsURI = att.getNamespaceURI();
            String attLocalName = att.getNodeName();

            // construct the FQN
            String attQName;
/*            if (attNsURI != null)
                attQName = attNsURI + attLocalName;
            else {
                attNsURI = attQName = "";
            }
*/
            if (attNsURI == null) attNsURI = "";
            attQName = attNsURI + attLocalName;

            String value = att.getNodeValue();

            attr.addAttribute(attNsURI, attLocalName, attQName, "", value); 
        }

        handler.startElement(nsURI, localName, qName, attr);
        handleChildren(node, handler);
        handler.endElement(nsURI, localName, qName);
    }

    static private void handleText(Node node, ContentHandler handler)
        throws SAXException
    {
        String text = node.getNodeValue();
        if (text != null) {
            char [] chars = text.toCharArray();
            handler.characters(chars, 0, text.length());
        }
    }

    static private void handleChildren(Node node, ContentHandler handler)
        throws SAXException
    {
        NodeList nodes = node.getChildNodes();

        for (int i=0; i < nodes.getLength(); i++) {
            convertToSAX(nodes.item(i), handler);
        }
    }

    public void setProperty(String prop, Object value)
        throws SAXNotRecognizedException
    {
        if (prop.equals(NAMESPACESURI)) {
            if (value instanceof Boolean) {
                supplyLocalName = ((Boolean) value).booleanValue();
            }
        }
        else if (prop.equals(NAMESPACEPREFIXURI)) {
            if (value instanceof Boolean) {
                supplyPrefix = ((Boolean) value).booleanValue();
            }
        }
        else throw new SAXNotRecognizedException("Property " + prop + " is not recognized.");
    }
}
