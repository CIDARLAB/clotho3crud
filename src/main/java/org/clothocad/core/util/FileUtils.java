/*
Copyright (c) 2009 The Regents of the University of California.
All rights reserved.
Permission is hereby granted, without written agreement and without
license or royalty fees, to use, copy, modify, and distribute this
software and its documentation for any purpose, provided that the above
copyright notice and the following two paragraphs appear in all copies
of this software.

IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY
FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES
ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
THE UNIVERSITY OF CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE
PROVIDED HEREUNDER IS ON AN "AS IS" BASIS, AND THE UNIVERSITY OF
CALIFORNIA HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES,
ENHANCEMENTS, OR MODIFICATIONS.
 */
package org.clothocad.core.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.Writer;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;


/**
 * @author sixpi
 */
public class FileUtils {
    final static Logger logger = LoggerFactory.getLogger(FileUtils.class);
    /**
     * dump an InputStream, returning a newly created byte[] array
     * @author Kelvin Li
     */
    public static byte[] dumpInputStream(InputStream instream) 
                                        throws IOException {
        int lastindex = 0;
        int increment;
        byte[] A = new byte[2048];
        byte[] B;
        do {
            B = new byte[2*A.length];
            System.arraycopy(A, 0, B, 0, A.length);
            A = B;
            try {
                increment = instream.read(A, lastindex, A.length - lastindex);
            } catch (IOException e) {
                logger.warn("cannot read InputStream", e);
                throw e;
            }
            if (increment == -1)
                /* EOF, no characters read */
                break;
            lastindex += increment;
        } while (lastindex >= A.length);

        /* copy A into B of precisely correct length */
        B = new byte[lastindex];
        System.arraycopy(A, 0, B, 0, B.length);
        return B;
    }

    /**
     * @param fileName
     * @return
     */
    public static boolean isXMLFile(String fileName) {
        return isFileType(fileName, "xml");
    }

    /**
     * @param fileName
     * @param extension without the .
     * @return
     */
    public static boolean isFileType(String fileName, String extension) {
        return fileName.endsWith( "." + extension );
    }

    public static String readFile(String path) {
        File f = new File(path);
        if (!f.isFile())
            return "";

        FileInputStream s;
        try {
            s = new FileInputStream(f);
        } catch (FileNotFoundException e) {
            logger.warn("", e);
            return "";
        }

        byte[] content;
        try {
            content = dumpInputStream(s);
        } catch (IOException e) {
            return "";
        }

        return new String(content);
    }

    public static void writeFile(String datafile, String filePath) {
        try {
            Writer output = null;
            File file = new File(filePath);
            output = new FileWriter(file);
            output.write(datafile);
            output.close();
        } catch (IOException ex) {
            ex.printStackTrace();
        }
    }

   /**
    * Convenience method to determine the filepath of the persisted data.
    * @param feature
    * @return 
    */
    public static String getFilePath(String fileName, String DBName) {
        File afile = new File(DBName);
        createDir(afile);
        return afile.getAbsolutePath() + File.separator + fileName;
    }

    /**
     * Helper method, creates directory with parent directories recursiverly.
     */
    private static void createDir(File f) {
        int limit = 10;
        while (!f.exists()) {
            if (!f.mkdir()) {
                createDir(f.getParentFile());
            }
            limit --;
            if(limit < 1) {
                break;
            }
        }
        if (limit == 0) {
            logger.warn("Directory was not created.");
        }
    }
}
