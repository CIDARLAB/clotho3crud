package org.clothocad.core.aspects.Interpreter;

import java.io.BufferedInputStream;
import java.util.Scanner;

public class StdIn {
    private String charsetName = "UTF-8";

    public Scanner scanner;

    public StdIn() {
        scanner = new Scanner(new BufferedInputStream(System.in), charsetName);
    }

    public String readString() {
        return scanner.nextLine();
    }

    public int readInt() {
        return scanner.nextInt();
    }
}
