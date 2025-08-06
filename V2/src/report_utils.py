import csv
from datetime import datetime

def generate_html_report(directory, output_html, combined_counts, summary_stats):
    """Generate HTML report"""
    html_path = os.path.join(directory, output_html)
    with open(html_path, 'w') as f:
        f.write(f"""
<!DOCTYPE html>
<html>
<head>
    <title>Read Recovery Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
        th, td {{ border: 1px solid #ddd; padding: 12px; text-align: left; }}
        th {{ background-color: #4CAF50; color: white; }}
        .number {{ font-family: monospace; }}
    </style>
</head>
<body>
    <h1>Read Recovery Report</h1>
    <p>Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
    
    <h2>Summary</h2>
    <table>
        <tr><th>Parameter</th><th>Value</th></tr>
        <tr><td>Max Mismatches</td><td class="number">{summary_stats['max_mismatches']}</td></tr>
        <tr><td>Reads Sampled</td><td class="number">{summary_stats['total_sampled']:,}</td></tr>
        <tr><td>Recoverable Reads</td><td class="number">{summary_stats['total_recoverable']:,}</td></tr>
        <tr><td>Recovery Rate</td><td class="number">{summary_stats['recovery_rate']:.2f}%</td></tr>
    </table>
    
    <h2>Sample Results</h2>
    <table>
        <tr><th>Sample</th><th>Original</th><th>Recovered</th><th>Total</th><th>Recovery %</th></tr>
        """)
        
        for sample, counts in sorted(combined_counts.items()):
            recovery_pct = (counts['recovered']/counts['total'])*100 if counts['total'] > 0 else 0
            f.write(f"""
        <tr>
            <td>{sample}</td>
            <td class="number">{counts['original']:,}</td>
            <td class="number">{counts['recovered']:,}</td>
            <td class="number">{counts['total']:,}</td>
            <td class="number">{recovery_pct:.2f}%</td>
        </tr>""")
        
        f.write("""
    </table>
</body>
</html>""")
    
    return html_path

def generate_csv_report(directory, output_csv, combined_counts):
    """Generate CSV report"""
    csv_path = os.path.join(directory, output_csv)
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['Sample', 'Original_Reads', 'Recovered_Reads', 'Total_Reads', 'Recovery_Percent'])
        for sample, counts in sorted(combined_counts.items()):
            recovery_pct = (counts['recovered']/counts['total'])*100 if counts['total'] > 0 else 0
            writer.writerow([sample, counts['original'], counts['recovered'], counts['total'], f"{recovery_pct:.2f}"])
    
    return csv_path