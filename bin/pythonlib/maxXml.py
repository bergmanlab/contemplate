#!/usr/bin/env python
import logging, re, urllib
import xml.etree.cElementTree as etree

class XmlParser:
    """ class to represent an xml tree (using ElementTree)
        Functions Accept PATH which is a /a/b/c style xpath-like expression to refer to elements
        PATH is not a complete XPATH implementation

        getText... functions return just a string
        getXml... functions return an XmlParser-object
        ...First  functions get only the first instance
        ...All    functions return an iterator
    
    >>> xp = XmlParser(string="<fruit><apple size='big'>boskoop</apple><apple size='small'>granny smith</apple><pear>mypear</pear></fruit>") 
    >>> xp.getTextFirst("pineapple", default="NothingAtAll")
    'NothingAtAll'
    >>> xp.getTextFirst("apple")
    'boskoop'
    >>> list(xp.getTextAll("apple"))
    ['boskoop', 'granny smith']
    >>> list(xp.getTextAll("apple", reqAttrDict={'size':'big'}))
    ['boskoop']
    
    """
    def __init__(self, string=None, url=None, root=None):
        self.root=None
        if string:
            self.fromString(string)
        elif url:
            self.fromUrl(url)
        elif root:
            self.root=root

    def getText(self):
        return self.root.text

    def fromString(self, string, removeNamespaces=False):
        root = etree.fromstring(string)
        if removeNamespaces:
            logging.debug("Stripping all namespaces")
            stripNamespaces(root)
        self.root = root

    def fromUrl(self, url, removeNamespaces=False, stopWords=[]):
        logging.debug("Retrieving %s" % url)
        text = urllib.urlopen(url).read()
        self.fromString(text, removeNamespaces=removeNamespaces)
        #for w in stopWords:
            #if w in text:
                #return None

    def _removeNamespaces(self):
        """ removes all namespaces from elementtree IN PLACE """
        root = self.root
        for el in root.getiterator():
            if el.tag[0] == '{':
                el.tag = el.tag.split('}', 1)[1]

    def _hasAttribs(self, el, reqAttrDict):
        for attr, value in reqAttrDict.iteritems():
            if el.attrib.get(attr, None)!=value:
                return False
        return True

    def getTextFirst(self, path, reqAttrDict=None, default=None):
        """ return text between elements given path 
            reqAttrDict is in the format attrName -> value
        """
        xml = self._getElFirst(path, reqAttrDict)
        if xml != None:
            return xml.text
        else:
            return default
        
    def getTextAll(self, path, reqAttrDict=None):
        for el in self._getElAll(path, reqAttrDict):
            yield el.text

    def _getElFirst(self, path, reqAttrDict):
        found = False
        for el in self._getElAll(path, reqAttrDict):
            found = True
            return el
        return None

    def _getElAll(self, path, reqAttrDict):
        found = False
        elIter = self.root.findall(path)
        for el in elIter:
            if reqAttrDict == None or self._hasAttribs(el, reqAttrDict):
                found = True
                yield el
        
    def getXmlFirst(self, path, reqAttrDict=None, default=None):
        el = self._getElFirst(path, reqAttrDict)
        if el==None:
            return default
        else:
            return XmlParser(root=el)

    def getXmlAll(self, path, reqAttrDict=None):
        for el in self._getElAll(path, reqAttrDict):
            yield XmlParser(root=el)

# ----- 
if __name__ == "__main__":
    import doctest
    doctest.testmod()
