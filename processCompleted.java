import java.awt.*;
import java.awt.event.*;
import javax.swing.*;
import java.awt.image.BufferedImage;
import java.io.File;
import javax.imageio.*;
import javax.imageio.stream.*;

public class processCompleted extends JFrame {

private static final long serialVersionUID = 1L;

public void ShowOK() {
       
try {
	
		final JFrame okFrame = new JFrame("");
        Container contentPane = okFrame.getContentPane();
        JButton okBtn = new JButton("OK");
        okBtn.setPreferredSize(new Dimension(30,30));
        okBtn.setBackground(Color.MAGENTA);
        okBtn.setForeground(Color.WHITE);
        
        String pathTofile; 
        pathTofile = System.getProperty("user.dir");
        System.out.print(pathTofile);
        pathTofile = pathTofile + "\\Omega.bmp";
        BufferedImage image;
        ImageInputStream fin = ImageIO.createImageInputStream( new File( pathTofile ) );
        image = ImageIO.read(new File( pathTofile ));
        fin.close();
        JLabel picture = new JLabel(new ImageIcon( image )) ;
        okFrame.add( picture );
        
        okFrame.setLocation(400,200);
        okFrame.setVisible(true);
        okFrame.setSize(100,100);
        contentPane.add(okBtn, BorderLayout.SOUTH);
        okFrame.pack();
        
        okBtn.addActionListener(new ActionListener() {       
       @Override public void actionPerformed(ActionEvent e) { 
        
            okFrame.dispose();
                     
        }});  
                 
        } catch (Exception x) {

          System.out.println(x);
          System.exit(0);
                              
        } 
}

public static void main(String[] args) {

(new processCompleted()).ShowOK();

}

}