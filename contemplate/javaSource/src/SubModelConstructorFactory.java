
import java.util.Map;
import java.util.HashMap;

import org.w3c.dom.NamedNodeMap;
import org.w3c.dom.Node;

public class SubModelConstructorFactory
{
    private static final Map constructorMap = new HashMap();

    // put code to register instances of ConstructorMaker here.
    static 
    {
        Block.registerMaker(constructorMap);
        Type1spacer.registerMaker(constructorMap);
        Type2spacer.registerMaker(constructorMap);
    }

    public static SubModelConstructor getConstructor(Node paramNode)
        throws IllegalArgumentException
    {
        String tag;
        NamedNodeMap attrs;

        // get the key
        if ((paramNode.getNodeType() == Node.ELEMENT_NODE) &&
            (tag = paramNode.getNodeName()).equals("submodel") &&
            ((attrs = paramNode.getAttributes()) != null)
            ) {
            String subModelKey;
            Node attr;
            SubModelConstructor.Maker maker;
            if (((attr = attrs.getNamedItem("type")) != null) &&
                ((subModelKey = attr.getNodeValue()) != null) &&
                ((maker = (SubModelConstructor.Maker) constructorMap.get(subModelKey)) != null)
                ) {
                return maker.makeConstructor(paramNode);
            }
        }
        throw new IllegalArgumentException("specified node was unsuccessfully parsed.");
    }
}

