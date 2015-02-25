
import java.util.Map;
import java.util.HashMap;

import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

public class ModelConstructorFactory
{
    private static final Map constructorMap = new HashMap();

    // put code to register instances of ConstructorMaker here.
    static 
    {
        EnhSearch.registerMaker(constructorMap);
    }

    public static ModelConstructor getConstructor(Node paramNode)
        throws IllegalArgumentException
    {
        String tag;
        NamedNodeMap attrs;

        // get the key
        if ((paramNode.getNodeType() == Node.ELEMENT_NODE) &&
            (tag = paramNode.getNodeName()).equals("markov_model") &&
            ((attrs = paramNode.getAttributes()) != null)
            ) {
            String subModelKey;
            Node attr;
            ModelConstructor.Maker maker;
            if (((attr = attrs.getNamedItem("type")) != null) &&
                ((subModelKey = attr.getNodeValue()) != null) &&
                ((maker = (ModelConstructor.Maker) constructorMap.get(subModelKey)) != null)
                ) {
                return maker.makeConstructor(paramNode);
            }
        }
        throw new IllegalArgumentException("specified node was unsuccessfully parsed.");
    }
}

