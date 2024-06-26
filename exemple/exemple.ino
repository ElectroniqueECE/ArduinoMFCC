#include <Arduino.h>
#include "arduinoMFCC.h"
// Définition des paramètres MFCC



#define MFCC_SIZE 11
#define DCT_MFCC_SIZE 6

#define FRAME_SIZE  256
#define FREQ_ECH 8000
// Déclaration de l'objet MFCC
arduinoMFCC mymfcc(MFCC_SIZE,DCT_MFCC_SIZE, FRAME_SIZE, FREQ_ECH);
//
//exemple de vecteur data
float frame[FRAME_SIZE]={
6.00,1591.00,1568.00,1708.00,1669.00,1680.00,1670.00,1614.00,1690.00,1641.00,1623.00,1588.00,1659.00,
1537.00,1382.00,1360.00,1304.00,1402.00,1221.00,1411.00,1447.00,1488.00,1506.00,1428.00,1467.00,1504.00
,1570.00,1659.00,1820.00,1891.00,1981.00,1933.00,1953.00,1846.00,1863.00,1778.00,1751.00,1758.00,1596.00
,1562.00,1628.00,1495.00,1409.00,1397.00,1412.00,1303.00,1479.00,1401.00,1608.00,1564.00,1614.00,1773.00,
1732.00,1838.00,1917.00,1822.00,1721.00,1730.00,1701.00,1574.00,1693.00,1543.00,1540.00,1397.00,1383.00,
1381.00,1373.00,1273.00,1474.00,1351.00,1487.00,1482.00,1517.00,1433.00,1637.00,1479.00,1524.00,1411.00,
1393.00,1456.00,1428.00,1437.00,1453.00,1518.00,1460.00,1524.00,1398.00,1372.00,1341.00,1395.00,1494.00,
1352.00,1366.00,1398.00,1526.00,1603.00,1632.00,1624.00,1678.00,1634.00,1766.00,1606.00,1673.00,1592.00,
1624.00,1537.00,1319.00,1345.00,1511.00,1624.00,1576.00,1618.00,1648.00,1511.00,1539.00,3.42};
// Déclaration de l'objet Audio
float mfcc[MFCC_SIZE];

void setup() {
  Serial.begin(9600);
  while(!Serial);

  mymfcc.create_hamming_window();
  mymfcc.create_mel_filter_bank(); 
  mymfcc.create_dct_matrix();
  


}

void loop() {
  
mymfcc.compute(frame, mfcc);

   for(int i=0 ;i<MFCC_SIZE;i++)
  {  Serial.println(mfcc[i]);

  }
  Serial.println("=======");

delay(1000);

  
}
